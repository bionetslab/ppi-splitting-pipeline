library(data.table)
library(ggplot2)
library(cowplot)
library(ggpattern)

# --precomputed_split, split seed 1234, sampling seed 42:
# candidate pool size (2x/4x/8x the positives) and the GO-BP weight (lambda_GO = 1 vs 0, at 4x).
repo <- "/path/to/project/ppi-splitting-pipeline/"
path_pool <- paste0(repo, "results_pool_size/multiqc/multiqc_report_data/")
path_go   <- paste0(repo, "results_go_weight/multiqc/multiqc_report_data/")

rename_map <- c(
  PDB_dimers = "PDB-Dimers",
  Pinder     = "PINDER",
  Hippie     = "HIPPIE",
  Intact     = "IntAct",
  String     = "STRING\nhuman, physical",
  String900  = "STRING-\n900",
  StringDB   = "STRING-\ndatabase",
  StringExp  = "STRING-\nexperimental",
  StringText = "STRING-\ntextmining"
)
rename_map_no_breaks <- c(
  PDB_dimers = "PDB-Dimers",
  Pinder     = "PINDER",
  Hippie     = "HIPPIE",
  Intact     = "IntAct",
  String     = "STRING: human, physical",
  String900  = "STRING-900",
  StringDB   = "STRING-database",
  StringExp  = "STRING-experimental",
  StringText = "STRING-textmining"
)
rename_split_map <- c(train = "Train", val = "Val.", test_balanced = "Test (bal.)")
rename_map_bias <- c(
  bias_deg_term  = "Degree",
  bias_self_term = "Self interactions",
  bias_tax_term  = "Taxonomy",
  bias_jac_term  = "GO-BP Jaccard index"
)
bias_colors <- c(
  "Degree"              = "#ffb000",
  "Self interactions"   = "#fe6100",
  "Taxonomy"            = "#dc267f",
  "GO-BP Jaccard index" = "#785ef0"
)

parse_ids <- function(dt) {
  dt[, `:=`(
    Dataset = sub("-.*", "", ID),
    Source  = factor(ifelse(grepl("-ILP-hcn-", ID), "HCN", "Random pool"), levels = c("Random pool", "HCN")),
    Config  = sub(".*-", "", ID)
  )]
}

read_diag <- function(path) {
  d <- fread(paste0(path, "multiqc_neg_ilp_diagnostics_table.txt"), sep = "\t")
  d <- melt(d, id.vars = c("ID", "split"), measure.vars = names(rename_map_bias), variable.name = "term")
  parse_ids(d)
  d[, `:=`(
    Dataset = factor(rename_map[Dataset], levels = unname(rename_map)),
    split   = factor(rename_split_map[split], levels = unname(rename_split_map)),
    term    = factor(rename_map_bias[as.character(term)], levels = unname(rename_map_bias))
  )]
}

#### Fig A3: objective terms vs candidate pool size ####

pool_diag <- read_diag(path_pool)
pool_diag[, pool_factor := as.integer(sub("pool", "", Config))]
# Taxonomy is weighted 0 for the single-species datasets, so drop terms that are 0 everywhere in a dataset
pool_diag <- pool_diag[, if (any(value != 0)) .SD, by = .(Dataset, term)]

# Two rows of datasets (5 + 4) instead of 9 columns; same y-range in both halves so the panels stay comparable
pool_panel <- function(ds) {
  ggplot(pool_diag[Dataset %in% ds], aes(x = pool_factor, y = value, color = term, linetype = Source, shape = Source)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.8) +
    scale_x_continuous(trans = "log2", breaks = c(2, 4, 8), labels = c("2×", "4×", "8×")) +
    scale_color_manual(values = bias_colors, drop = FALSE) +
    coord_cartesian(ylim = range(pool_diag$value)) +
    facet_grid(split ~ Dataset) +
    labs(x = "Candidate pool size (× positives)", y = "Objective value", color = NULL, linetype = NULL, shape = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          strip.text.x = element_text(face = "bold", size = 10, margin = margin(b = 4)),
          strip.text.y = element_text(face = "bold", size = 10))
}
ds_levels <- levels(pool_diag$Dataset)
p_top    <- pool_panel(ds_levels[1:5]) + theme(axis.title.x = element_blank())
p_bottom <- pool_panel(ds_levels[6:9])
pool_legend <- cowplot::get_legend(pool_panel(ds_levels) + theme(legend.position = "bottom"))
plot_grid(p_top, plot_grid(p_bottom, NULL, rel_widths = c(4.2, 0.8)), pool_legend,
          ncol = 1, rel_heights = c(1, 1, 0.08))
ggsave("figures/FigureA3.pdf", width = 10, height = 10, dpi = 300)

#### Fig A4: Random Forest AUPRC ####

config_labels <- c(pool2 = "2× pool", pool4 = "4× pool (default)", go0 = "4× pool, no GO-BP term", pool8 = "8× pool")
config_colors <- c("2× pool" = "#ffb000", "4× pool (default)" = "#785ef0", "8× pool" = "#dc267f",
                   "4× pool, no GO-BP term" = "#648fff")

read_rf <- function(path) {
  rbindlist(lapply(c("balanced", "realistic"), function(ts)
    fread(paste0(path, "multiqc_classifier_metrics_test_", ts, "_table.txt"), sep = "\t")[, TestSet := ts]))
}
rf <- rbind(read_rf(path_pool), read_rf(path_go))
parse_ids(rf)
rf[, `:=`(
  Dataset = factor(rename_map_no_breaks[Dataset], levels = rev(unname(rename_map_no_breaks))),
  Config  = factor(config_labels[Config], levels = rev(unname(config_labels)))
)]
dodge_width <- position_dodge(width = 0.9)
rf_panel <- function(ts, baseline, xlim, hjust, reverse) {
  p <- ggplot(rf[TestSet == ts], aes(y = Dataset, x = AUPRC, fill = Config, pattern = Source)) +
    geom_col_pattern(position = dodge_width, pattern_fill = "black", pattern_colour = NA,
                     pattern_density = 0.25, pattern_spacing = 0.015, pattern_angle = 45) +
    geom_text(aes(label = sprintf("%.2f", AUPRC)), position = dodge_width, hjust = hjust, size = 2.6) +
    geom_vline(xintercept = baseline, linetype = 2) +
    facet_grid(Source ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = config_colors, breaks = unname(config_labels)) +
    scale_pattern_manual(values = c("Random pool" = "none", "HCN" = "stripe")) +
    guides(fill = guide_legend(override.aes = list(pattern = "none"), nrow = 2),
           pattern = guide_legend(override.aes = list(fill = "grey80"))) +
    theme_minimal() +
    theme(legend.position = "none", axis.title.y = element_blank()) +
    ggtitle(if (ts == "balanced") "Balanced" else "Realistic")
  if (reverse) {
    p + scale_x_reverse() + coord_cartesian(xlim = xlim, expand = FALSE) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.text.y = element_blank(),
            plot.margin = margin(5.5, 2, 5.5, 5.5))
  } else {
    p + coord_cartesian(xlim = xlim, expand = FALSE) + theme(plot.margin = margin(5.5, 5.5, 5.5, 2))
  }
}

p_balanced  <- rf_panel("balanced", 0.5, c(1.0, 0.45), 1.05, reverse = TRUE)
p_realistic <- rf_panel("realistic", 1/11, c(0.0, 1.0), -0.1, reverse = FALSE)

shared_legend <- cowplot::get_legend(p_realistic + labs(fill = NULL, pattern = NULL) + theme(legend.position = "bottom"))
plots_row <- plot_grid(p_balanced, p_realistic, nrow = 1, align = "h", axis = "tb")
final_plot <- plot_grid(plots_row, shared_legend, ncol = 1, rel_heights = c(1, 0.08))
final_plot
ggsave("figures/FigureA4.pdf", final_plot, width = 9, height = 8, dpi = 300)

#### Fig A5: bias NMI (as in Fig 4a) with vs without the GO-BP term (4x pool) ####

read_bias <- function(root) {
  f <- Sys.glob(paste0(repo, root, "/*/multiqc/*_bias_mqc.tsv"))
  d <- rbindlist(lapply(f, fread), idcol = "i")
  d[, `:=`(ID = basename(dirname(dirname(f[i]))), Bias = sub("_bias_mqc.tsv", "", basename(f[i])), i = NULL)]
}
# test_realistic negatives come from the degree sampler, not the ILP, so that split is identical in both runs
nmi <- rbind(read_bias("results_pool_size"), read_bias("results_go_weight"))[Bias != "topology_shortcut" & Split != "test_realistic"]
parse_ids(nmi)
nmi <- nmi[Config %in% c("pool4", "go0")]
nmi_split_map <- c(train = "Train", val = "Validation", test_balanced = "Test (balanced)")
nmi_bias_map <- c(self_interactions = "Self interactions", same_species = "Same Species",
                  sequence_similarity = "Sequence Similarity", embedding_similarity = "Embedding Similarity",
                  functional_relatedness_BP = "GO-BP", functional_relatedness_MF = "GO-MF",
                  functional_relatedness_CC = "GO-CC")
nmi[, `:=`(
  Dataset = factor(rename_map[Dataset], levels = unname(rename_map)),
  Split   = factor(nmi_split_map[Split], levels = unname(nmi_split_map)),
  Bias    = factor(nmi_bias_map[Bias], levels = unname(nmi_bias_map)),
  Config  = factor(c(pool4 = "With GO-BP term", go0 = "Without GO-BP term")[Config],
                   levels = c("With GO-BP term", "Without GO-BP term"))
)]

ggplot(nmi, aes(x = Config, y = `NMI(A;Y)`, color = Bias, linetype = Source, shape = Source,
                group = interaction(Bias, Source))) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Self interactions" = "#ffb000", "Sequence Similarity" = "#fe6100",
                                "Same Species" = "#fc92d7", "Embedding Similarity" = "#c701ff",
                                "GO-BP" = "#648fff", "GO-MF" = "#003a7d", "GO-CC" = "#bdd373")) +
  facet_grid(Split ~ Dataset) +
  guides(color = guide_legend(nrow = 1)) +
  labs(y = "NMI(a; y)", x = NULL, color = NULL, linetype = NULL, shape = NULL) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "top",
        strip.text.x = element_text(face = "bold", size = 12, margin = margin(b = 4)),
        strip.text.y = element_text(face = "bold", size = 12, margin = margin(b = 4)),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/FigureA5.pdf", width = 15, height = 8, dpi = 300)
