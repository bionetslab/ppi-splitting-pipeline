library(data.table)
library(ggplot2)
library(dplyr)

core_path <- "/path/to/ppi-splitting-pipeline/results_prott5/"
core_path_others <- "/path/to/ppi-splitting-pipeline/results_server/"
path_to_prott5_results <- paste0(core_path, "multiqc_prott5/")
path_to_results <- paste0(core_path_others, "multiqc_non_struc/")
path_to_struc_results <- paste0(core_path_others, "multiqc_structural/")
path_to_pdb_results <- paste0(core_path_others, "multiqc_PDB/")
path_to_hcn_results <- paste0(core_path_others, "multiqc_hcn/")


##### Bias analysis ######
rename_map <- c(
  Hippie     = "HIPPIE",
  Intact     = "IntAct",
  String     = "STRING\nhuman, physical",
  String900  = "STRING-\n900",
  StringDB  = "STRING-\ndatabase",
  StringExp  = "STRING-\nexperimental",
  StringText  = "STRING-\ntextmining",
  PDB_dimers  = "PDB-Dimers",
  Pinder  = "PINDER"
)

all_biases <- lapply(list.files(core_path, pattern="*mqc.tsv", full.names=TRUE, recursive=TRUE), fread)
names(all_biases) <- list.files(core_path, pattern="*mqc.tsv", full.names=TRUE, recursive=TRUE)

all_biases_esm2 <- lapply(list.files(core_path_others, pattern="*mqc.tsv", full.names=TRUE, recursive=TRUE), fread)
names(all_biases_esm2) <- list.files(core_path_others, pattern="*mqc.tsv", full.names=TRUE, recursive=TRUE)

all_biases <- rbindlist(all_biases, idcol="filename")
all_biases_esm2 <- rbindlist(all_biases_esm2, idcol="filename")
all_biases$embedding <- "ProtT5"
all_biases_esm2$embedding <- "ESM-2"

all_biases[, filename := gsub(paste0(core_path, "/"), "", filename, fixed=TRUE)]
all_biases_esm2[, filename := gsub(paste0(core_path_others, "/"), "", filename, fixed=TRUE)]
all_biases <- rbind(all_biases, all_biases_esm2)

all_biases[, filename := gsub("multiqc/", "", filename, fixed=TRUE)]
all_biases[, filename := gsub("_bias_mqc.tsv", "", filename, fixed=TRUE)]
all_biases[, c("Sample", "Bias") := tstrsplit(filename, "/")]
all_biases[, filename := NULL]
all_biases$Dataset <- sapply(stringr::str_split(all_biases$Sample, "-", n = 2),`[`, 1)
all_biases$Category <- sapply(stringr::str_split(all_biases$Sample, "-", n = 2),`[`, 2)

only_leakage <- all_biases[Category == "Leakage-split"]
only_leakage[Dataset %in% names(rename_map), Dataset := rename_map[Dataset]]
only_leakage[, Dataset := factor(Dataset, levels = c("PDB-Dimers", 
                                                     "PINDER", 
                                                     "IntAct", 
                                                     "HIPPIE", 
                                                     "STRING\nhuman, physical", 
                                                     "STRING-\n900",
                                                     "STRING-\ndatabase",
                                                     "STRING-\nexperimental",
                                                     "STRING-\ntextmining"))]

rename_split_map <- c(
  train = "Train",
  val = "Validation",
  test_balanced = "Test (balanced)",
  test_realistic = "Test (realistic)"
)
only_leakage[Split %in% names(rename_split_map), Split := rename_split_map[Split]]
only_leakage[, Split := factor(Split, levels = c("Train", "Validation", "Test (balanced)", "Test (realistic)"))]

rename_bias_map <- c(
  self_interactions = "Self interactions",
  same_species = "Same Species",
  topology_shortcut = "Topology",
  sequence_similarity = "Sequence Similarity",
  embedding_similarity = "Embedding Similarity",
  functional_relatedness_BP = "GO-BP",
  functional_relatedness_MF = "GO-MF",
  functional_relatedness_CC = "GO-CC"
)
only_leakage[Bias %in% names(rename_bias_map), Bias := rename_bias_map[Bias]]
only_leakage[, Bias := factor(Bias, levels = c("Self interactions", 
                                               "Same Species", 
                                               "Topology",
                                               "Sequence Similarity",
                                               "Embedding Similarity",
                                               "GO-BP",
                                               "GO-MF",
                                               "GO-CC"))]
only_leakage[`Detectability (Spearman ρ)` == "NaN", `Detectability (Spearman ρ)` := NA]
only_leakage <- only_leakage[!(Bias == "Same Species" & !Dataset %in% c("PDB-Dimers", "PINDER", "IntAct"))]

comparison <- dcast(only_leakage, ... ~ embedding, value.var = c("NMI(A;Y)", "Detectability (Spearman ρ)"))
comparison[,`NMI ESM2 - NMI ProtT5` := `NMI(A;Y)_ESM-2` - `NMI(A;Y)_ProtT5`]
comparison[,`Detectability ESM2 - Detectability ProtT5` := `Detectability (Spearman ρ)_ESM-2` - `Detectability (Spearman ρ)_ProtT5`]

ggplot(only_leakage[Bias == "Embedding Similarity"], aes(
  x = embedding, 
  y = `NMI(A;Y)`, color = Split, group = Split))+
  geom_point()+
  geom_line()+
  geom_hline(aes(yintercept = 0.0))+
  scale_color_manual(values = c(
    "Train" = "#ffb000",
    "Validation" = "#fe6100",
    "Test (balanced)" = "#dc267f",
    "Test (realistic)" = "#785ef0"
  )) +
  facet_wrap(~Dataset)+
  labs(x = "Embedding", y = "Utility")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90))
ggsave("figures/prott5_nmi_comp_leakage.pdf", width=9, height=4)

ggplot(comparison, aes(x = Bias, y = `Detectability ESM2 - Detectability ProtT5`, color = Split))+
  geom_point()+
  geom_hline(aes(yintercept = 0.0))+
  facet_wrap(~Dataset, ncol=5)+
  scale_color_manual(values = c(
    "Train" = "#ffb000",
    "Validation" = "#fe6100",
    "Test (balanced)" = "#dc267f",
    "Test (realistic)" = "#785ef0"
  )) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90))
ggsave("figures/prott5_detectability_comp_leakage.pdf", width=8, height=5)

### Random Forest comparison ####
test_results_balanced_prott5 <- fread(paste0(path_to_prott5_results, "multiqc_report_data/multiqc_classifier_metrics_test_balanced_table.txt"), sep="\t")
test_results_balanced_prott5$embedding <- "ProtT5"

test_results_balanced <- fread(paste0(path_to_results, "multiqc_report_data/multiqc_classifier_metrics_test_balanced_table.txt"), sep="\t")
test_results_balanced_struc <- fread(paste0(path_to_struc_results, "multiqc_report_data/multiqc_classifier_metrics_test_balanced_table.txt"), sep="\t")
test_results_balanced_struc <- test_results_balanced_struc[!startsWith(Sample, "PDB")]
test_results_balanced_pdb <- fread(paste0(path_to_pdb_results, "multiqc_report_data/multiqc_classifier_metrics_test_balanced_table.txt"), sep="\t")
test_results_balanced_hcn <- fread(paste0(path_to_hcn_results, "multiqc_report_data/multiqc_classifier_metrics_test_balanced_table.txt"), sep="\t")
test_results_balanced <- rbind(test_results_balanced, test_results_balanced_struc)
test_results_balanced <- rbind(test_results_balanced, test_results_balanced_pdb)
test_results_balanced$embedding <- "ESM-2"

test_results_balanced <- rbind(test_results_balanced, test_results_balanced_prott5)
test_results_balanced$Dataset <- sapply(stringr::str_split(test_results_balanced$Sample, "-", n = 2),`[`, 1)
test_results_balanced$Category <- sapply(stringr::str_split(test_results_balanced$Sample, "-", n = 2),`[`, 2)
test_results_balanced[, Category := factor(Category, levels = c("Leakage-split", "regular", "ILP-split", "ILP-split-and-sampling", "ILP-hcn"))]

test_results_real_prott5 <- fread(paste0(path_to_prott5_results, "multiqc_report_data/multiqc_classifier_metrics_test_realistic_table.txt"), sep="\t")
test_results_real_prott5$embedding <- "ProtT5"

test_results_real <- fread(paste0(path_to_results, "multiqc_report_data/multiqc_classifier_metrics_test_realistic_table.txt"), sep="\t")
test_results_real_struc <- fread(paste0(path_to_struc_results, "multiqc_report_data/multiqc_classifier_metrics_test_realistic_table.txt"), sep="\t")
test_results_real_struc <- test_results_real_struc[!startsWith(Sample, "PDB")]
test_results_real_pdb <- fread(paste0(path_to_pdb_results, "multiqc_report_data/multiqc_classifier_metrics_test_realistic_table.txt"), sep="\t")
test_results_real_hcn <- fread(paste0(path_to_hcn_results, "multiqc_report_data/multiqc_classifier_metrics_test_realistic_table.txt"), sep="\t")
test_results_real <- rbind(test_results_real, test_results_real_struc)
test_results_real <- rbind(test_results_real, test_results_real_pdb)
test_results_real <- rbind(test_results_real, test_results_real_hcn)
test_results_real$embedding <- "ESM-2"

test_results_real <- rbind(test_results_real, test_results_real_prott5)
test_results_real$Dataset <- sapply(stringr::str_split(test_results_real$Sample, "-", n = 2),`[`, 1)
test_results_real$Category <- sapply(stringr::str_split(test_results_real$Sample, "-", n = 2),`[`, 2)
test_results_real[, Category := factor(Category, levels = c("Leakage-split", "regular", "ILP-split", "ILP-split-and-sampling", "ILP-hcn"))]

test_results_balanced$TestSet <- "balanced"
test_results_real$TestSet <- "realistic"

rename_map <- c(
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

comparison <- rbind(test_results_balanced, test_results_real)
comparison <- comparison[!(Dataset == "Intact" & Category == "ILP-hcn")]
comparison[Dataset %in% names(rename_map), Dataset := rename_map[Dataset]]
comparison[, Dataset := factor(Dataset, levels = rev(c("PDB-Dimers", 
                                                       "PINDER", 
                                                       "IntAct", 
                                                       "HIPPIE", 
                                                       "STRING: human, physical", 
                                                       "STRING-900",
                                                       "STRING-experimental",
                                                       "STRING-database",
                                                       "STRING-textmining")))]

comparison <- comparison[, -c("Sample", "ID")]

category_map <- c(
  "Leakage-split" = "Maximally biased",
  "regular" = "Similarity-reduced",
  "ILP-split" = "ILP split",
  "ILP-split-and-sampling" = "ILP split + ILP sampling",
  "ILP-hcn" = "ILP split + ILP sampling (HCN)"
)

comparison[Category %in% names(category_map), Category := category_map[Category]]
comparison[, Category := factor(Category, levels = rev(c(
  "Maximally biased", "Similarity-reduced", "ILP split", "ILP split + ILP sampling", "ILP split + ILP sampling (HCN)"
)))]

comparison <- comparison[, c("AUPRC", "embedding", "Dataset", "Category", "TestSet")]
comparison <- dcast(comparison, ... ~ embedding, value.var = "AUPRC")
comparison[, AUPRC_diff := `ESM-2` - `ProtT5`]

dodge_width <- position_dodge(width = 0.9)
library(cowplot)
p_balanced <- ggplot(comparison %>% filter(TestSet == "balanced"),
                     aes(y = Dataset, x = AUPRC_diff, fill = Category)) +
  geom_col(position = dodge_width) +
  geom_text(aes(label = paste0(round(`ESM-2`, 2), "-", round(`ProtT5`, 2)),
                hjust = ifelse(AUPRC_diff >= 0, -0.1, 1.1)),
            position = dodge_width,
            size = 3.2) +
  scale_fill_manual(values = c(
    "Maximally biased" = "#ffb000",
    "Similarity-reduced" = "#fe6100",
    "ILP split" = "#dc267f",
    "ILP split + ILP sampling" = "#785ef0", 
    "ILP split + ILP sampling (HCN)" = "#648fff"
  )) +                                 
  coord_cartesian(xlim = c(-0.4, 0.4), expand = FALSE) +
  theme_minimal() +
  labs(x = "AUPRC ESM-2 - AUPRC ProtT5")+
  theme(legend.position = "none",
        axis.text.y  = element_blank(),   
        axis.ticks.y = element_blank(),   
        axis.title.y = element_blank(),   
        plot.margin = margin(5.5, 2, 5.5, 5.5)) +
  ggtitle("Balanced")

p_realistic <- ggplot(comparison %>% filter(TestSet == "realistic"), 
                      aes(y = Dataset, x = AUPRC_diff, fill = Category)) +
  geom_col(position = dodge_width) +
  geom_text(aes(label = paste0(round(`ESM-2`, 2), "-", round(`ProtT5`, 2)),
                hjust = ifelse(AUPRC_diff >= 0, -0.1, 1.1)),
            position = dodge_width,
            size = 3.2) +
  scale_fill_manual(values = c(
    "Maximally biased" = "#ffb000",
    "Similarity-reduced" = "#fe6100",
    "ILP split" = "#dc267f",
    "ILP split + ILP sampling" = "#785ef0",
    "ILP split + ILP sampling (HCN)" = "#648fff"
  )) +
  coord_cartesian(xlim = c(-0.4, 0.4), expand = FALSE) +
  theme_minimal() +
  labs(x = "AUPRC ESM-2 - AUPRC ProtT5")+
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
ggsave("figures/prott5_rf_comparison.pdf", height=6, width=8)

