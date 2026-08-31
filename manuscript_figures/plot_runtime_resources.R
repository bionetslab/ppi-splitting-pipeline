library(data.table)
library(ggplot2)
library(stringr)
library(cowplot)
library(scales)

report_paths <- c(
  prott5_run = "/path/to/report-20260818-57891661.html"
)

run_labels <- c(
  prott5_run = "ProtT5 Run"
)

parse_nf_report <- function(path) {
  content <- paste(readLines(path, warn = FALSE), collapse = "\n")

  data_start <- str_locate(content, fixed("window.data = {"))[1, "end"]
  body <- str_sub(content, data_start + 1)

  data_end <- str_locate(body, "\\]\\s*\\}\\s*;")[1, "end"]
  body <- str_sub(body, 1, data_end)

  rec_starts <- gregexpr('\\{"task_id"', body, perl = TRUE)[[1]]
  rec_ends <- c(rec_starts[-1] - 1, nchar(body))
  records <- substring(body, rec_starts, rec_ends)

  get_field <- function(key) str_match(records, paste0('"', key, '":"([^"]*)"'))[, 2]

  data.table(
    process  = get_field("process"),
    tag      = get_field("tag"),
    status   = get_field("status"),
    realtime = as.numeric(get_field("realtime")),
    peak_rss = as.numeric(get_field("peak_rss")),
    cpus     = as.numeric(get_field("cpus"))
  )
}

all_tasks <- rbindlist(lapply(names(report_paths), function(nm) {
  dt <- parse_nf_report(path.expand(report_paths[[nm]]))
  dt[, run := run_labels[[nm]]]
  dt
}))

all_tasks[, `:=`(
  module   = sub(":.*$", "", process),
  runtime_min = realtime / 1000 / 60,
  peak_mem_gb = peak_rss / 1024^3
)]

### Readable labels #############################################################

module_labels <- c(
  DATA_PREP        = "Data preparation",
  CLUSTERING       = "Clustering",
  SPLIT_POSITIVES  = "Positive split",
  SAMPLE_NEGATIVES = "Negative sampling",
  BIAS_DIAGNOSTICS = "Quality control",
  QC               = "Quality control",
  TRAIN_BASELINE   = "Baseline model"
)

process_labels <- c(
  "DATA_PREP:FETCH_DATA"                     = "Fetch Data (UniProt)",
  "DATA_PREP:GET_LENGTHS_SHARED"             = "Compute sequence lengths",
  "DATA_PREP:GET_LENGTHS"                    = "Compute sequence lengths",
  "DATA_PREP:SUBSET_FETCHED_DATA"            = "Subset fetched data",
  "CLUSTERING:RUN_BLAST"                     = "All-vs-all BLAST",
  "CLUSTERING:MAKE_METIS"                    = "Build METIS graph",
  "CLUSTERING:RUN_KAHIP"                     = "KaHIP clustering",
  "SPLIT_POSITIVES:SOLVE_ILP"                = "Solve ILP (split)",
  "SPLIT_POSITIVES:CDHIT2D"                  = "CD-HIT2D redundancy filter",
  "SPLIT_POSITIVES:REMOVE_REDUNDANT"         = "Remove redundant pairs",
  "SPLIT_POSITIVES:SORT_PPIS"                = "Sort PPIs",
  "SPLIT_POSITIVES:SPLIT_RANDOM"             = "Random split",
  "SAMPLE_NEGATIVES:SAMPLE_NEGATIVES_ILP"    = "Solve ILP (negative sampling)",
  "SAMPLE_NEGATIVES:SAMPLE_NEGATIVES_DEGREE" = "Degree-based sampling",
  "BIAS_DIAGNOSTICS:BIAS_ANALYSIS"           = "Bias analysis",
  "QC:COLLECT_BIAS"                          = "Collect bias metrics",
  "QC:SIMILARITY_HEATMAP"                    = "Similarity heatmap",
  "QC:MULTIQC"                               = "MultiQC report",
  "TRAIN_BASELINE:EMBED_SEQUENCES"           = "Compute embeddings",
  "TRAIN_BASELINE:TRAIN_CLASSIFIER"          = "Train baseline classifier"
)

all_tasks[, module_label := module_labels[module]]
all_tasks[, process_label := process_labels[process]]

# KaHIP is run with a different number of target clusters depending on the
# splitting strategy (k=3 for regular/leakage splits, k=100 for ILP-based
# ones), which strongly affects its runtime, so keep the two apart.
all_tasks[, kahip_k := str_match(tag, ":\\s*k=(\\d+)$")[, 2]]
all_tasks[process == "CLUSTERING:RUN_KAHIP",
          process_label := paste0("KaHIP clustering (k=", kahip_k, ")")]

order_tbl <- all_tasks[, .(max = max(runtime_min)), by = process_label][order(max)]
all_tasks[, process_label := factor(process_label, levels = order_tbl$process_label)]
all_tasks[, module_label := factor(module_label, levels = c("Data preparation", "Clustering", 
                                                            "Positive split", "Negative sampling", 
                                                            "Baseline model", "Quality control"))]

module_colors <- c(
  "Data preparation"   = "#bababa",
  "Clustering"          = "#fe6100",
  "Positive split"       = "#dc267f",
  "Negative sampling"    = "#785ef0",
  "Quality control"      = "#648fff",
  "Baseline model"       = "#ffb000"
)

### Plot #########################################################################

p_time <- ggplot(all_tasks, aes(x = process_label, y = runtime_min, color = module_label)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.3) +
  #scale_y_log10(labels = label_number()) +
  scale_color_manual(values = module_colors) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(x = NULL, y = "Wall-clock runtime [min]", color = NULL) +
  theme(legend.position = "none")

p_mem <- ggplot(all_tasks, aes(x = process_label, y = peak_mem_gb, color = module_label)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.3) +
  #scale_y_log10(labels = label_number()) +
  scale_color_manual(values = module_colors) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(x = NULL, y = "Peak memory [GB]", color = NULL) +
  theme(axis.text.y = element_blank())

legend_source <- p_mem + theme(legend.position = "top")+
  guides(color = guide_legend(nrow = 1))
shared_legend <- cowplot::get_legend(legend_source)

plots_row <- plot_grid(p_time, p_mem + theme(legend.position = "none"),
                        nrow = 1, rel_widths = c(1.35, 1))
final_plot <- plot_grid(shared_legend, plots_row, ncol = 1, rel_heights = c(0.08, 1))
final_plot

ggsave("figures/runtime_resources_prott5.pdf", final_plot, width = 8, height = 5, dpi = 300)

### Summary table for the manuscript text ######################################

run_summary <- all_tasks[, .(
  n_tasks             = .N,
  total_task_hours    = sum(realtime) / 1000 / 3600,
  total_cpu_hours     = sum(realtime * cpus, na.rm = TRUE) / 1000 / 3600,
  max_peak_mem_gb     = max(peak_mem_gb)
), by = run]
print(run_summary)

step_summary <- all_tasks[, .(
  median_min = median(runtime_min),
  max_min    = max(runtime_min),
  median_gb  = median(peak_mem_gb),
  max_gb     = max(peak_mem_gb)
), by = .(module_label, process_label)][order(-median_min)]
print(step_summary)

fwrite(step_summary, "runtime_resource_summary.csv")
