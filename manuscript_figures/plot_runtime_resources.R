library(data.table)
library(ggplot2)
library(stringr)
library(cowplot)
library(scales)

repo <- "/path/to/ppi-splitting-pipeline/"
# The ESM-2 robustness runs, concatenated from several Nextflow reports (missing ones are skipped)
report_paths <- paste0(repo, "results_robustness/", c("runtime_robustness.html", "runtime_robustness_struc.html",
                                                     "runtime_robustness_hcn.html", "runtime_robustness_bias.html",
                                                     "runtime_robustness_hcn_bias.html",
                                                     "robustness_rerun.html", "robustness_bias_rerun.html"))

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

if (any(!file.exists(report_paths))) message("Skipping missing reports: ", paste(basename(report_paths[!file.exists(report_paths)]), collapse = ", "))
all_tasks <- rbindlist(lapply(report_paths[file.exists(report_paths)], function(p) parse_nf_report(p)[, report := basename(p)]))
# The ILP-split-and-sampling rows were rerun (--precomputed_split) after the negative-sampling fix. Take their tasks
# from negative sampling on from the rerun reports; clustering/splitting wasn't rerun, so those stay from the old ones.
# The rerun's other tasks (EXTRACT_POSITIVES, GET_LENGTHS, MULTIQC) only exist because of --precomputed_split.
# A missing rerun report keeps the old tasks.
rerun_reports <- c("robustness_rerun.html", "robustness_bias_rerun.html")
rerun_procs <- intersect(c("SAMPLE_NEGATIVES:SAMPLE_NEGATIVES_ILP",
                           "SAMPLE_NEGATIVES:SAMPLE_NEGATIVES_DEGREE", "TRAIN_BASELINE:TRAIN_CLASSIFIER",
                           "ROBUSTNESS:ROBUSTNESS_SUMMARY", "BIAS_DIAGNOSTICS:BIAS_ANALYSIS"),
                         all_tasks[report %in% rerun_reports, process])
all_tasks <- all_tasks[ifelse(report %in% rerun_reports, process %in% rerun_procs,
                              !(process %in% rerun_procs & grepl("-ILP-split-and-sampling", tag)))]
# Failed attempts are retried; only count the attempts that produced results (cached tasks keep their original metrics)
all_tasks <- all_tasks[status != "FAILED"]

# The ILP negative sampling in the HCN report ran with cvxpy 1.9.2, whose Gurobi interface builds the model
# in O(n^2) Python (hours for the large train splits). Replace it with the default-configuration runs (4x pool, all
# lambdas as in the robustness runs) from the --precomputed_split pool-size run with cvxpy >= 1.9.3, on the released
# splits (split seed 1234, sampling seed 42).
# The random-pool rows already ran with the fixed cvxpy in the rerun report, so those stay when present.
fixed_ilp_report <- paste0(repo, "results_pool_size/pool_size.html")
if (file.exists(fixed_ilp_report)) {
  fixed_ilp <- parse_nf_report(fixed_ilp_report)[process == "SAMPLE_NEGATIVES:SAMPLE_NEGATIVES_ILP" & status != "FAILED"]
  fixed_ilp[, report := basename(fixed_ilp_report)]
  has_rerun_ilp <- "SAMPLE_NEGATIVES:SAMPLE_NEGATIVES_ILP" %in% rerun_procs
  all_tasks <- rbind(
    all_tasks[process != "SAMPLE_NEGATIVES:SAMPLE_NEGATIVES_ILP" | report %in% rerun_reports],
    fixed_ilp[grepl(if (has_rerun_ilp) "-ILP-hcn-pool4_" else "-ILP-(hcn-)?pool4_", tag)]
  )
} else {
  message("Missing ", basename(fixed_ilp_report), ": ILP negative sampling still shows the slow cvxpy 1.9.2 runtimes")
}

# Steps the ESM-2 robustness runs don't have, from the ProtT5 run: the embedding step (ESM-2 embeddings were
# precomputed; one task embeds every protein) and the normal-mode QC subworkflow (robustness mode skips it)
prott5_report <- paste0(repo, "results_prott5/report-prott5.html")
all_tasks <- rbind(all_tasks, parse_nf_report(prott5_report)[
  (process == "TRAIN_BASELINE:EMBED_SEQUENCES" | grepl("^QC:", process)) & status != "FAILED"][
  , report := basename(prott5_report)])

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
  ROBUSTNESS       = "Quality control",
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
  "ROBUSTNESS:ROBUSTNESS_SUMMARY"            = "Robustness summary",
  "ROBUSTNESS:MULTIQC"                       = "MultiQC report (robustness)",
  "TRAIN_BASELINE:EMBED_SEQUENCES"           = "Compute embeddings (ProtT5)",
  "TRAIN_BASELINE:TRAIN_CLASSIFIER"          = "Train baseline classifier"
)

all_tasks[, module_label := module_labels[module]]
all_tasks[, process_label := process_labels[process]]

# KaHIP is run with a different number of target clusters depending on the
# splitting strategy (k=3 for regular/leakage splits, k=100 for ILP-based
# ones), which strongly affects its runtime, so keep the two apart.
all_tasks[, kahip_k := str_match(tag, ":\\s*k=(\\d+)")[, 2]]  # tag ends in ", seed=<split seed>" in newer runs
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

### Plot ########################################################################

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
print(final_plot)

ggsave("figures/FigureC1.pdf", final_plot, width = 8, height = 5, dpi = 300)

### Summary table for the manuscript text ######################################

# Some cached tasks lack peak memory, hence na.rm for the memory columns
run_summary <- all_tasks[, .(
  n_tasks             = .N,
  total_task_hours    = sum(realtime) / 1000 / 3600,
  total_cpu_hours     = sum(realtime * cpus, na.rm = TRUE) / 1000 / 3600,
  max_peak_mem_gb     = max(peak_mem_gb, na.rm = TRUE)
)]
print(run_summary)

step_summary <- all_tasks[, .(
  median_min = median(runtime_min),
  max_min    = max(runtime_min),
  median_gb  = median(peak_mem_gb, na.rm = TRUE),
  max_gb     = max(peak_mem_gb, na.rm = TRUE)
), by = .(module_label, process_label)][order(-median_min)]
print(step_summary)

fwrite(step_summary, "runtime_resource_summary.csv")
