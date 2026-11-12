library(data.table)
library(ggplot2)
library(dplyr)

# Figs B1/B2: ProtT5 bias-only run on the same Leakage-split robustness splits as the ESM-2 one (1234/42), with ProtT5
# embeddings restricted to proteins that have an ESM-2 embedding (samplesheet_prott5_same_splits_bias_only.csv)
core_path <- "/path/to/ppi-splitting-pipeline/results_prott5_same_splits_bias/"
core_path_others <- "/path/to/ppi-splitting-pipeline/results_robustness_bias/"  # ESM-2 bias-only run on the 1234/42 robustness splits
# Fig B3: ProtT5 RF retrained on the ESM-2 robustness splits (split 1234, sampling 42, training 42) on the same machine,
# so ESM-2 and ProtT5 differ only in the embeddings (bin/retrain_on_esm2_splits.py)
path_to_prott5_results <- "/path/to/ppi-splitting-pipeline/results_prott5_same_splits/"


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
ggsave("figures/FigureB1.pdf", width=9, height=4)

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
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
ggsave("figures/FigureB2.pdf", width=8, height=5)

### Random Forest comparison ####
test_results_balanced_prott5 <- fread(paste0(path_to_prott5_results, "multiqc_classifier_metrics_test_balanced_table.txt"), sep="\t")
test_results_balanced_prott5$embedding <- "ProtT5"

# ESM-2: robustness runs with split seed 1234, sampling seed 42, training seed 42 (the pipeline defaults)
path_to_robustness <- "/path/to/ppi-splitting-pipeline/results_robustness/"
path_to_hcn <- "/path/to/ppi-splitting-pipeline/results_robustness_hcn/"
esm2_runs <- rbindlist(lapply(Sys.glob(paste0(c(path_to_robustness, path_to_hcn), "*/robustness/robustness_runs.tsv")), fread, sep = "\t"))
esm2_runs <- esm2_runs[split_seed == 1234 & sampling_seed == 42 & training_seed == 42]
esm2_runs[, `:=`(Sample = ID, embedding = "ESM-2")]

test_results_balanced <- rbind(esm2_runs[test_set == "test_balanced", names(test_results_balanced_prott5), with = FALSE],
                               test_results_balanced_prott5)
test_results_balanced$Dataset <- sapply(stringr::str_split(test_results_balanced$Sample, "-", n = 2),`[`, 1)
test_results_balanced$Category <- sapply(stringr::str_split(test_results_balanced$Sample, "-", n = 2),`[`, 2)
test_results_balanced[, Category := factor(Category, levels = c("Leakage-split", "regular", "ILP-split", "ILP-split-and-sampling", "ILP-hcn"))]

test_results_real_prott5 <- fread(paste0(path_to_prott5_results, "multiqc_classifier_metrics_test_realistic_table.txt"), sep="\t")
test_results_real_prott5$embedding <- "ProtT5"

test_results_real <- rbind(esm2_runs[test_set == "test_realistic", names(test_results_real_prott5), with = FALSE],
                           test_results_real_prott5)
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
  "ILP-hcn" = "ILP split + ILP sampling (HRNI)"
)

comparison[Category %in% names(category_map), Category := category_map[Category]]
comparison[, Category := factor(Category, levels = rev(c(
  "Maximally biased", "Similarity-reduced", "ILP split", "ILP split + ILP sampling", "ILP split + ILP sampling (HRNI)"
)))]

# Table B1: ProtT5 metrics, laid out like rf_table.tex (plot_comparison.R); one seed combination, so no SD
library(kableExtra)
library(stringr)
metrics <- c("AUROC", "AUPRC", "F1", "MCC", "Precision", "Recall", "Accuracy")
prott5_dt <- comparison[embedding == "ProtT5", c("Dataset", "Category", "TestSet", metrics), with = FALSE][
  order(-Dataset, -Category, TestSet)][, (metrics) := lapply(.SD, function(x) sprintf("%.2f", x)), .SDcols = metrics]
prott5_dt[, Category := as.character(Category)]
prott5_dt[Category == "Maximally biased", Category := "Maximally\nbiased"]
prott5_dt[Category == "Similarity-reduced", Category := "Similarity-\nreduced"]
prott5_dt[Category == "ILP split + ILP sampling (HRNI)", Category := "ILP split + sampling (HRNI)"]  # shorter, fits the column

# Dataset becomes a group-header row (\multicolumn{9}{l}{Dataset: \textbf{...}}) instead of a column
dataset_runs <- rle(as.character(prott5_dt$Dataset))
dataset_index <- setNames(dataset_runs$lengths, paste0("Dataset: \\textbf{", dataset_runs$values, "}"))

prott5_dt %>%
  mutate(Category = linebreak(str_replace_all(Category, " \\+ ", " +\n"), align = "tl")) %>%  # "tl" = top + left
  select(Category, TestSet, all_of(metrics)) %>%
  kbl(format = "latex", booktabs = TRUE, longtable = TRUE,
      linesep = "", escape = FALSE,
      col.names = c("Category", "Test set", metrics),
      caption = "Performance metrics of the baseline with ProtT5 embeddings across datasets, splitting strategies, and test sets.") %>%
  kable_styling(latex_options = c("repeat_header")) %>%
  collapse_rows(columns = 1, latex_hline = "none", valign = "top") %>%  # must run before pack_rows, else it mangles \textbf
  column_spec(1, width = "1.8cm") %>%  # after collapse_rows: the other order drops a } from \makecell cells
  pack_rows(index = dataset_index, escape = FALSE, bold = FALSE, indent = FALSE, latex_gap_space = "0em") %>%
  # \midrule above and below each Dataset row; pack_rows only offers \hline, which clashes with booktabs
  gsub("\\\\addlinespace\\[0em\\]\n", "\\\\midrule\n", .) %>%
  gsub("(\\\\multicolumn\\{[0-9]+\\}\\{l\\}\\{Dataset:[^\n]*\n)", "\\1\\\\midrule\n", .) %>%
  sub("\\\\endlastfoot\n\\\\midrule\n", "\\\\endlastfoot\n", .) %>%  # first group already sits under the header rule
  paste0("{\\small\\setlength{\\tabcolsep}{3pt}\n", ., "}") %>%  # braces keep font/spacing local to this table
  writeLines("figures/rf_table_prott5.tex")

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
    "ILP split + ILP sampling (HRNI)" = "#648fff"
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
    "ILP split + ILP sampling (HRNI)" = "#648fff"
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
ggsave("figures/FigureB3.pdf", height=7, width=9)


### Fig 5a (plot_hcn.R) with ProtT5: ILP split + sampling vs. HRNI ####
# Same retrained runs as Fig B3; one seed combination, so no error bars
hrni_dataset_map <- rename_map[c("Hippie", "String", "String900", "StringDB", "StringExp", "StringText")]
hrni_category_map <- c(
  "ILP-split-and-sampling" = "ILP-based split and sampling",
  "ILP-hcn" = "ILP-based split and sampling from HRNI"
)
hrni <- rbindlist(list(balanced = test_results_balanced_prott5, realistic = test_results_real_prott5), idcol = "TestSet")
# IntAct is multi-species, the HRNI network is human-only; PDB-Dimers/PINDER have no HRNI run
hrni <- hrni[grepl("-(ILP-split-and-sampling|ILP-hcn)$", ID) & sub("-.*", "", ID) %in% names(hrni_dataset_map)]
hrni[, `:=`(
  Dataset  = factor(hrni_dataset_map[sub("-.*", "", ID)], levels = rev(unname(hrni_dataset_map))),
  Category = factor(hrni_category_map[sub("^[^-]*-", "", ID)], levels = rev(unname(hrni_category_map)))
)]

hrni_fill <- scale_fill_manual(values = c(
  "ILP-based split and sampling" = "#785ef0",
  "ILP-based split and sampling from HRNI" = "#648fff"
))

p5_balanced <- ggplot(hrni[TestSet == "balanced"], aes(y = Dataset, x = AUPRC, fill = Category)) +
  geom_col(position = dodge_width) +
  geom_text(aes(label = sprintf("%.2f", AUPRC)), position = dodge_width, hjust = 1.05, size = 3.2) +
  geom_vline(xintercept = 0.5, linetype = 2) +
  hrni_fill +
  scale_x_reverse() +
  coord_cartesian(xlim = c(1.0, 0.48), expand = FALSE) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(5.5, 2, 5.5, 5.5)) +
  ggtitle("Balanced")

p5_realistic <- ggplot(hrni[TestSet == "realistic"], aes(y = Dataset, x = AUPRC, fill = Category)) +
  geom_col(position = dodge_width) +
  geom_text(aes(label = sprintf("%.2f", AUPRC)), position = dodge_width, hjust = -0.1, size = 3.2) +
  geom_vline(xintercept = 1/11, linetype = 2) +
  hrni_fill +
  coord_cartesian(xlim = c(0.0, 1.0), expand = FALSE) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = margin(5.5, 5.5, 5.5, 2)) +
  ggtitle("Realistic")

p5_legend <- cowplot::get_legend(p5_realistic + labs(fill = NULL) + theme(legend.position = "bottom"))
p5_row <- plot_grid(p5_balanced, p5_realistic, nrow = 1, align = "h", axis = "tb", rel_widths = c(1, 1))
final_plot <- plot_grid(p5_row, p5_legend, ncol = 1, rel_heights = c(1, 0.1))
final_plot
ggsave("figures/fig5a_hcn_rf_prott5.pdf", final_plot, width = 8, height = 3, dpi = 300)
