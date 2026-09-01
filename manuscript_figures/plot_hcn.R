library(data.table)
library(ggplot2)
library(cowplot)
library(dplyr)

### Fig. 5a: Random Forest results ####
core_path <- "/path/to/ppi-splitting-pipeline/results_esm2/"

path_to_results <- paste0(core_path, "multiqc_non_struc/")
path_to_hcn_results <- paste0(core_path, "multiqc_hcn/")

test_results_balanced <- fread(paste0(path_to_results, "multiqc_report_data/multiqc_classifier_metrics_test_balanced_table.txt"), sep="\t")
test_results_balanced_hcn <- fread(paste0(path_to_hcn_results, "multiqc_report_data/multiqc_classifier_metrics_test_balanced_table.txt"), sep="\t")
test_results_balanced <- rbind(test_results_balanced, test_results_balanced_hcn)
test_results_balanced$Dataset <- sapply(stringr::str_split(test_results_balanced$Sample, "-", n = 2),`[`, 1)
test_results_balanced$Category <- sapply(stringr::str_split(test_results_balanced$Sample, "-", n = 2),`[`, 2)
test_results_balanced[, Category := factor(Category, levels = c("Leakage-split", "regular", "ILP-split", "ILP-split-and-sampling", "ILP-hcn"))]

test_results_real <- fread(paste0(path_to_results, "multiqc_report_data/multiqc_classifier_metrics_test_realistic_table.txt"), sep="\t")
test_results_real_hcn <- fread(paste0(path_to_hcn_results, "multiqc_report_data/multiqc_classifier_metrics_test_realistic_table.txt"), sep="\t")
test_results_real <- rbind(test_results_real, test_results_real_hcn)
test_results_real$Dataset <- sapply(stringr::str_split(test_results_real$Sample, "-", n = 2),`[`, 1)
test_results_real$Category <- sapply(stringr::str_split(test_results_real$Sample, "-", n = 2),`[`, 2)
test_results_real[, Category := factor(Category, levels = c("Leakage-split", "regular", "ILP-split", "ILP-split-and-sampling", "ILP-hcn"))]

test_results_balanced$TestSet <- "balanced"
test_results_real$TestSet <- "realistic"

comparison <- rbind(test_results_balanced, test_results_real)
rename_map_no_breaks <- c(
  Hippie     = "HIPPIE",
  Intact     = "IntAct",
  String     = "STRING: human, physical",
  String900  = "STRING-900",
  StringDB  = "STRING-database",
  StringExp  = "STRING-experimental",
  StringText  = "STRING-textmining",
  PDB_dimers  = "PDB-Dimers",
  Pinder  = "PINDER"
)
comparison <- comparison[!Dataset %in% c("Intact", "PDB_dimers", "Pinder"), ]
comparison[, Dataset := rename_map_no_breaks[Dataset]]
comparison[, Dataset := factor(Dataset, levels = rev(c("HIPPIE", 
                                                       "STRING: human, physical", 
                                                       "STRING-900",
                                                       "STRING-database",
                                                       "STRING-experimental",
                                                       "STRING-textmining")))]

comparison <- comparison[, -c("Sample", "ID")]

category_map <- c(
  "Leakage-split" = "Maximally biased",
  "regular" = "Similarity-reduced",
  "ILP-split" = "ILP-based split",
  "ILP-split-and-sampling" = "ILP-based split and sampling",
  "ILP-hcn" = "ILP-based split and sampling from HCN"
)

comparison[Category %in% names(category_map), Category := category_map[Category]]
comparison[, Category := factor(Category, levels = rev(c(
  "Maximally biased", "Similarity-reduced", "ILP-based split", "ILP-based split and sampling",  "ILP-based split and sampling from HCN"
)))]

comparison <- comparison[Category %in% c("ILP-based split and sampling",  "ILP-based split and sampling from HCN")]

dodge_width <- position_dodge(width = 0.9)

p_balanced <- ggplot(comparison %>% filter(TestSet == "balanced"), 
                     aes(y = Dataset, x = AUPRC, fill = Category)) +
  geom_col(position = dodge_width) +
  geom_text(aes(label = sprintf("%.2f", AUPRC)),
            position = dodge_width,
            hjust = 1,       
            size = 3.2) +
  geom_vline(xintercept = 0.5, linetype = 2) +
  scale_fill_manual(values = c(
    "ILP-based split and sampling" = "#785ef0",
    "ILP-based split and sampling from HCN" = "#648fff"
  )) +
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
  geom_text(aes(label = sprintf("%.2f", AUPRC)),
            position = dodge_width,
            hjust = -0.1,
            size = 3.2) +
  geom_vline(xintercept = 0.1, linetype = 2) +
  scale_fill_manual(values = c(
    "ILP-based split and sampling" = "#785ef0",
    "ILP-based split and sampling from HCN" = "#648fff"
  )) +
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

residuals <- fread(paste0(path_to_results, "multiqc_report_data/multiqc_neg_ilp_diagnostics_table.txt"), sep="\t")
residuals <- residuals[, c("ID", "split", "obj_value", "bias_deg_term", "bias_tax_term", "bias_self_term", "bias_jac_term", "status")]
res_hcn <- fread(paste0(path_to_hcn_results, "multiqc_report_data/multiqc_neg_ilp_diagnostics_table.txt"), sep="\t")
res_hcn <- res_hcn[, c("ID", "split", "obj_value", "bias_deg_term", "bias_tax_term", "bias_self_term", "bias_jac_term", "status")]
residuals <- rbind(residuals, res_hcn)
residuals$Dataset <- sapply(stringr::str_split(residuals$ID, "-", n = 2),`[`, 1)
residuals$Category <- sapply(stringr::str_split(residuals$ID, "-", n = 2),`[`, 2)
rename_map_cat <- c(
  "ILP-split-and-sampling" = "ILP",
  "ILP-hcn" = "ILP + HCN"
)
residuals[, Category := rename_map_cat[Category]]
residuals[, Category := factor(Category, levels = c("ILP", "ILP + HCN"))]
residuals <- residuals[Dataset != "Intact"]

res_long <- melt(residuals, id.vars = c("Dataset", "Category", "ID", "split", "obj_value", "status"))
res_long[split %in% names(rename_split_map), split := rename_split_map[split]]
res_long[, split := factor(split, levels = c("Train", "Validation", "Test (balanced)", "Test (realistic)"))]
res_long[Dataset %in% names(rename_map), Dataset := rename_map[Dataset]]
res_long[, Dataset := factor(Dataset, levels = c("PDB-Dimers", 
                                                 "PINDER", 
                                                 "IntAct", 
                                                 "HIPPIE", 
                                                 "STRING\nhuman, physical", 
                                                 "STRING-\n900",
                                                 "STRING-\ndatabase",
                                                 "STRING-\nexperimental",
                                                 "STRING-\ntextmining"))]
res_long <- res_long[variable != "bias_tax_term"]
rename_map_bias <- c(
  bias_deg_term = "Degree",
  bias_tax_term = "Taxonomy",
  bias_self_term = "Self interactions",
  bias_jac_term = "GO-BP Jaccard index"
)
res_long[, variable := rename_map_bias[variable]]
res_long[, variable := factor(variable, levels = c(
  "Degree",
  "Self interactions",
  "GO-BP Jaccard index"
))]

ggplot(res_long, aes(x = Category, y = value, fill = variable))+
  geom_col()+
  facet_grid(split~Dataset)+
  scale_fill_manual(values = c(
    "Degree"              = "#ffb000",
    "Self interactions"             = "#fe6100",
    #"Taxonomy"                 = "#dc267f",
    "GO-BP Jaccard index" = "#785ef0"
  )) +
  theme_minimal(base_size = 12)+
  labs(fill = NULL, y = "Objective value", x = NULL)+
  theme(legend.position = "bottom",
        strip.text.x  = element_text(face = "bold", size = 12, margin = margin(b = 4)),
        )
ggsave("figures/fig5b_objective_values.pdf", width=8, height=5, dpi=300)
