library(data.table)
library(ggplot2)
library(cowplot)
library(dplyr)

path_to_robustness <- "/path/to/project/ppi-splitting-pipeline/results_robustness/"
path_to_hcn <- "/path/to/project/ppi-splitting-pipeline/results_robustness_hcn/"

pop_sd <- function(x) sqrt(mean((x - mean(x))^2))  # population SD, as in bin/robustness_summary.py

# IntAct is multi-species, the HCN is human-only; PDB-Dimers/PINDER have no HCN run
keep_ids <- function(ids) grepl("-(ILP-split-and-sampling|ILP-hcn)$", ids) &
  !sub("-.*", "", ids) %in% c("Intact", "PDB_dimers", "Pinder")

category_map <- c(
  "ILP-split-and-sampling" = "ILP-based split and sampling",
  "ILP-hcn" = "ILP-based split and sampling from HCN"
)

### Fig. 5a: Random Forest results ####
rename_map_no_breaks <- c(
  Hippie     = "HIPPIE",
  String     = "STRING: human, physical",
  String900  = "STRING-900",
  StringDB  = "STRING-database",
  StringExp  = "STRING-experimental",
  StringText  = "STRING-textmining"
)

# Seed robustness runs (split x sampling x training seeds)
robustness_runs <- rbindlist(lapply(Sys.glob(paste0(c(path_to_robustness, path_to_hcn), "*/robustness/robustness_runs.tsv")),
                                    fread, sep = "\t"))
robustness_runs <- robustness_runs[keep_ids(ID)]
robustness_runs[, `:=`(
  Dataset  = factor(rename_map_no_breaks[sub("-.*", "", ID)], levels = rev(unname(rename_map_no_breaks))),
  Category = factor(category_map[sub("^[^-]*-", "", ID)], levels = rev(unname(category_map))),
  TestSet  = sub("test_", "", test_set)
)]

# Bars show the mean over seed runs
comparison <- robustness_runs[, .(
  AUPRC      = mean(AUPRC),
  AUPRC_sd   = pop_sd(AUPRC),
  AUPRC_max  = max(AUPRC)
), by = .(Dataset, Category, TestSet)]
comparison[, label_x := pmax(AUPRC + AUPRC_sd, AUPRC_max)]  # labels beyond error bars

dodge_width <- position_dodge(width = 0.9)
hcn_fill <- scale_fill_manual(values = c(
  "ILP-based split and sampling" = "#785ef0",
  "ILP-based split and sampling from HCN" = "#648fff"
))
error_bars <- geom_errorbar(aes(xmin = AUPRC - AUPRC_sd, xmax = AUPRC + AUPRC_sd),
                            position = dodge_width, width = 0.5, linewidth = 0.4, show.legend = FALSE)

p_balanced <- ggplot(comparison %>% filter(TestSet == "balanced"),
                     aes(y = Dataset, x = AUPRC, fill = Category)) +
  geom_col(position = dodge_width) +
  error_bars +
  geom_text(aes(x = label_x, label = sprintf("%.2f ± %.2f", AUPRC, AUPRC_sd)),
            position = dodge_width,
            hjust = 1.05,
            size = 3.2) +
  geom_vline(xintercept = 0.5, linetype = 2) +
  hcn_fill +
  scale_x_reverse() +
  coord_cartesian(xlim = c(1.0, 0.48), expand = FALSE)+
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(5.5, 2, 5.5, 5.5)
        )+
  ggtitle("Balanced")

p_realistic <- ggplot(comparison %>% filter(TestSet == "realistic"),
                      aes(y = Dataset, x = AUPRC, fill = Category)) +
  geom_col(position = dodge_width) +
  error_bars +
  geom_text(aes(x = label_x, label = sprintf("%.2f ± %.2f", AUPRC, AUPRC_sd)),
            position = dodge_width,
            hjust = -0.1,
            size = 3.2) +
  geom_vline(xintercept = 1/11, linetype = 2) +
  hcn_fill +
  coord_cartesian(xlim = c(0.0, 1.0), expand = FALSE) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = margin(5.5, 5.5, 5.5, 2))+
  ggtitle("Realistic")

# Build a throwaway plot just to extract its legend
legend_source <- p_realistic + labs(fill = NULL)+
  theme(legend.position="bottom")
shared_legend <- cowplot::get_legend(legend_source)

# Combine the two panels side by side
plots_row <- plot_grid(p_balanced, p_realistic, nrow = 1,
                       align = "h", axis = "tb",
                       rel_widths = c(1, 1))

# Stack the row of plots on top of the shared legend
final_plot <- plot_grid(plots_row, shared_legend, ncol = 1, rel_heights = c(1, 0.1))
final_plot
ggsave("figures/fig5a_hcn_rf.pdf", width=8, height=3, dpi=300)

### Fig 5b: Bias objective terms results ###

rename_map <- c(
  Hippie     = "HIPPIE",
  String     = "STRING\nhuman, physical",
  String900  = "STRING-\n900",
  StringDB  = "STRING-\ndatabase",
  StringExp  = "STRING-\nexperimental",
  StringText  = "STRING-\ntextmining"
)
rename_split_map <- c(
  train = "Train",
  val = "Val.",
  test_balanced = "Test (bal.)"
)

# Robustness runs, all 3x3 (split, sampling) seed combinations; these tables come from bin/reconstruct_neg_ilp_diagnostics.py
ilp_runs <- rbindlist(lapply(c(path_to_robustness, path_to_hcn), function(d)
  fread(paste0(d, "multiqc/multiqc_report_data/multiqc_neg_ilp_diagnostics_table.txt"), sep="\t")))
ilp_runs <- ilp_runs[keep_ids(ID)]
# Bars show the mean of each term over seed combinations, error bars its SD
# Taxonomy is left out: all remaining datasets are human-only
res_long <- melt(ilp_runs, id.vars = c("ID", "split"),
                 measure.vars = c("bias_deg_term", "bias_self_term", "bias_jac_term"))
res_long <- res_long[, .(value = mean(value), sd = pop_sd(value)), by = .(ID, split, variable)]
res_long[, `:=`(Dataset  = factor(rename_map[sub("-.*", "", ID)], levels = unname(rename_map)),
                Category = factor(c("ILP-split-and-sampling" = "ILP", "ILP-hcn" = "ILP + HCN")[sub("^[^-]*-", "", ID)],
                                  levels = c("ILP", "ILP + HCN")),
                split    = factor(rename_split_map[split], levels = unname(rename_split_map)))]
rename_map_bias <- c(
  bias_deg_term = "Degree",
  bias_self_term = "Self interactions",
  bias_jac_term = "GO-BP Jaccard index"
)
res_long[, variable := factor(rename_map_bias[variable], levels = unname(rename_map_bias))]

ggplot(res_long, aes(x = Category, y = value, fill = variable))+
  geom_col(position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9),
                width = 0.3, linewidth = 0.4)+
  facet_grid(split~Dataset)+
  scale_fill_manual(values = c(
    "Degree"              = "#ffb000",
    "Self interactions"             = "#fe6100",
    "GO-BP Jaccard index" = "#785ef0"
  )) +
  theme_minimal(base_size = 12)+
  labs(fill = NULL, y = "Objective value", x = NULL)+
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.text.x  = element_text(face = "bold", size = 12, margin = margin(b = 4)),
        )
ggsave("figures/fig5b_objective_values.pdf", width=8, height=5, dpi=300)
