library(data.table)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(fmsb)
library(dplyr)
library(scales)
library(cowplot)
library(gridGraphics)
library(Matrix)

core_path <- "/path/to/ppi-splitting-pipeline/results_esm2/"
path_to_results <- paste0(core_path, "multiqc_non_struc/")
path_to_struc_results <- paste0(core_path, "multiqc_structural/")
path_to_pdb_results <- paste0(core_path, "multiqc_PDB/")

##### Fig 2a: Bias analysis ######
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

all_biases <- rbindlist(all_biases, idcol="filename")
all_biases[, filename := gsub(paste0(core_path, "/"), "", filename, fixed=TRUE)]
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


ggplot(only_leakage, aes(y=`NMI(A;Y)`, x=`Detectability (Spearman ρ)`, shape=Split, color=Bias, group=Bias))+
  geom_point(size=2.5)+
  geom_path(linewidth = 0.4, alpha = 0.2) +
  scale_color_manual(values=c("Self interactions" = "#ffb000", 
                              "Same Species"="#fe6100", 
                              "Topology"="#fc92d7",
                              "Sequence Similarity" = "#c701ff",
                              "Embedding Similarity" = "#648fff",
                              "GO-BP" = "#003a7d",
                              "GO-MF" = "#d83034",
                              "GO-CC" = "#bdd373"))+
  scale_x_continuous(breaks = c(0, 0.5, 1.0)) +
  labs(color = NULL, shape = NULL) +
  facet_wrap(~Dataset, ncol=2)+
  theme_minimal(base_size = 16)+
  guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        legend.spacing.x = unit(0.1, "cm"),
        legend.margin = margin(0, 0, -10, 0) )


ggplot(only_leakage, aes(y=`NMI(A;Y)`, x=`Detectability (Spearman ρ)`, shape=Split, color=Dataset, group=Dataset))+
  geom_point(size=2.5)+
  geom_path(linewidth = 0.4, alpha = 0.2) +
  scale_x_continuous(breaks = c(0, 0.5, 1.0)) +
  scale_color_manual(values=c("PDB-Dimers" = "#ffb000", 
                              "PINDER"="#fe6100", 
                              "IntAct"="#fc92d7",
                              "HIPPIE" = "#c701ff",
                              "STRING\nhuman, physical" = "#648fff",
                              "STRING-\n900" = "#003a7d",
                              "STRING-\ndatabase" = "#d83034",
                              "STRING-\nexperimental" = "#bdd373",
                              "STRING-\ntextmining" = "#bababa"))+
  labs(color = NULL, shape = NULL) +
  facet_wrap(~Bias, ncol=4)+
  theme_minimal(base_size = 16)+
  labs(x = "Detectability: Ridge Regressor", y = "Utility: NMI(a; y)")+
  guides(color = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
  theme(panel.grid.minor = element_blank(), legend.position = "right",
        legend.spacing.x = unit(0.1, "cm"),
        legend.margin = margin(0, 0, 0, 0) )
ggsave("figures/fig2a_bias_analysis_leakage.pdf", width = 12, height = 5, dpi = 300)


#### Figure 2b: Radar performance ####

test_results_balanced <- fread(paste0(path_to_results, "multiqc_report_data/multiqc_classifier_metrics_test_balanced_table.txt"), sep="\t")
test_results_balanced_struc <- fread(paste0(path_to_struc_results, "multiqc_report_data/multiqc_classifier_metrics_test_balanced_table.txt"), sep="\t")
test_results_balanced_struc <- test_results_balanced_struc[!startsWith(Sample, "PDB")]
test_results_balanced_pdb <- fread(paste0(path_to_pdb_results, "multiqc_report_data/multiqc_classifier_metrics_test_balanced_table.txt"), sep="\t")
test_results_balanced <- rbind(test_results_balanced, test_results_balanced_struc)
test_results_balanced <- rbind(test_results_balanced, test_results_balanced_pdb)
test_results_balanced$Dataset <- sapply(stringr::str_split(test_results_balanced$Sample, "-", n = 2),`[`, 1)
test_results_balanced$Category <- sapply(stringr::str_split(test_results_balanced$Sample, "-", n = 2),`[`, 2)
test_results_balanced[, Category := factor(Category, levels = c("Leakage-split", "regular", "ILP-split", "ILP-split-and-sampling"))]

test_results_real <- fread(paste0(path_to_results, "multiqc_report_data/multiqc_classifier_metrics_test_realistic_table.txt"), sep="\t")
test_results_real_struc <- fread(paste0(path_to_struc_results, "multiqc_report_data/multiqc_classifier_metrics_test_realistic_table.txt"), sep="\t")
test_results_real_struc <- test_results_real_struc[!startsWith(Sample, "PDB")]
test_results_real_pdb <- fread(paste0(path_to_pdb_results, "multiqc_report_data/multiqc_classifier_metrics_test_realistic_table.txt"), sep="\t")
test_results_real <- rbind(test_results_real, test_results_real_struc)
test_results_real <- rbind(test_results_real, test_results_real_pdb)
test_results_real$Dataset <- sapply(stringr::str_split(test_results_real$Sample, "-", n = 2),`[`, 1)
test_results_real$Category <- sapply(stringr::str_split(test_results_real$Sample, "-", n = 2),`[`, 2)
test_results_real[, Category := factor(Category, levels = c("Leakage-split", "regular", "ILP-split", "ILP-split-and-sampling"))]

test_results_balanced$TestSet <- "balanced"
test_results_real$TestSet <- "realistic"

only_leakage <- rbind(test_results_balanced, test_results_real)
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
only_leakage[Dataset %in% names(rename_map_no_breaks), Dataset := rename_map_no_breaks[Dataset]]
only_leakage[, Dataset := factor(Dataset, levels = c("PDB-Dimers", 
                                                     "PINDER", 
                                                     "IntAct", 
                                                     "HIPPIE", 
                                                     "STRING: human, physical", 
                                                     "STRING-900",
                                                     "STRING-database",
                                                     "STRING-experimental",
                                                     "STRING-textmining"))]

only_leakage <- only_leakage[Category == "Leakage-split"]
only_leakage <- only_leakage[, -c("Sample", "ID", "Category")]
only_leakage <- only_leakage[order(AUROC)]
only_leakage <- only_leakage[order(Dataset, TestSet)]

metrics <- c("AUROC", "AUPRC", "F1", "MCC", "Precision", "Recall", "Accuracy")
n <- length(metrics)

get_radar_coords <- function(values, n) {
  angles <- seq(90, 90+360, length.out = n + 1)[1:n] * pi / 180
  data.frame(x = values * cos(angles), y = values * sin(angles))
}

make_panel <- function(ds) {
  par(mar = c(0.5, 0.5, 2.5, 0.5))
  
  sub  <- only_leakage %>% filter(Dataset == ds)
  bal  <- sub %>% filter(TestSet == "balanced")
  real <- sub %>% filter(TestSet == "realistic")
  radar_df <- rbind(
    # maximum values AUROC, AUPRC, F1, MCC, Precision, Recall, Accuracy
    rep(1, n),
    # maximum values AUROC, AUPRC, F1, MCC, Precision, Recall, Accuracy
    rep(0, n),
    as.numeric(bal[, ..metrics]),
    as.numeric(real[, ..metrics])
  )
  colnames(radar_df) <- metrics
  radar_df <- as.data.frame(radar_df)
  
  radarchart(radar_df,
             axistype = 0,
             pcol = c("#785ef0", "#fe6100"),
             pfcol = c(alpha("#785ef0", 0.15), alpha("#fe6100", 0.15)),
             plwd = 2.5, plty = 1,
             cglcol = "grey80", cglty = 1, cglwd = 0.5,
             vlcex = 1.5,                    
             title = ds,
             cex.main = 1.6)                  
  
  bal_coords  <- get_radar_coords(as.numeric(bal[, ..metrics]), n)
  real_coords <- get_radar_coords(as.numeric(real[, ..metrics]), n)
  
  text(bal_coords$x, bal_coords$y,
       labels = sprintf("%.2f", as.numeric(bal[, ..metrics])),
       col = "#785ef0", cex = 0.95, pos = 1)  
  
  text(real_coords$x, real_coords$y,
       labels = sprintf("%.2f", as.numeric(real[, ..metrics])),
       col = "#fe6100", cex = 0.95, pos = 3) 
}

datasets <- unique(only_leakage$Dataset)

panels <- lapply(datasets, function(ds) {
  eval(bquote(~ make_panel(.(ds))))
})

legend_panel <- ~ {
  plot.new()
  legend("center", legend = c("balanced", "realistic"),
         col = c("#785ef0", "#fe6100"), lwd = 2, bty = "n", cex = 1.8)
}

final_plot <- plot_grid(plotlist = c(panels, list(legend_panel)),
                        nrow = 2,
                        rel_widths = rep(1, 5),
                        rel_heights = rep(1, 2),
                        scale = 1)

ggsave("figures/fig2b_radar.pdf", final_plot, width = 20, height = 8, dpi = 300)

#### Figure 2c: Degree ####

leakage_dirs <- c("Hippie", "Intact", "PDB_dimers", "Pinder", "String", "String900",
                  "StringDB", "StringExp", "StringText")

split_files <- c(
  "Train" = "train.csv",
  "Validation" = "val.csv",
  "Test (balanced)" = "test_balanced.csv",
  "Test (realistic)" = "test_realistic.csv"
)

load_leakage_split <- function(ds) {
  base <- paste0(core_path, ds, "-Leakage-split/")
  dt <- rbindlist(lapply(names(split_files), function(split_name) {
    d <- fread(paste0(base, split_files[[split_name]]), select = c("protein1", "protein2", "label"))
    d[, Split := split_name]
    d
  }))
  dt[, Dataset := ds]
  dt
}

all_interactions <- rbindlist(lapply(leakage_dirs, load_leakage_split))

all_interactions[, rowid := .I]
protein_long <- melt(all_interactions, id.vars = c("rowid", "Dataset", "Split", "label"),
                     measure.vars = c("protein1", "protein2"), value.name = "protein")
# a self-interaction (protein1 == protein2) melts into two identical rows; dedup so it counts once
protein_long <- unique(protein_long, by = c("rowid", "protein"))

degree_dt <- protein_long[, .N, by = c("Dataset", "Split", "label", "protein")]
setnames(degree_dt, "N", "degree")

degree_dt[, label_lab := factor(ifelse(label == 1, "Positive", "Negative"), levels = c("Positive", "Negative"))]
degree_dt[, Split := factor(Split, levels = c("Train", "Validation", "Test (balanced)", "Test (realistic)"))]

degree_dt[Dataset %in% names(rename_map), Dataset := rename_map[Dataset]]
degree_dt[, Dataset := factor(Dataset, levels = c("PDB-Dimers",
                                                  "PINDER",
                                                  "IntAct",
                                                  "HIPPIE",
                                                  "STRING\nhuman, physical",
                                                  "STRING-\n900",
                                                  "STRING-\ndatabase",
                                                  "STRING-\nexperimental",
                                                  "STRING-\ntextmining"))]

ggplot(degree_dt[degree <=100 & Split == "Train"], aes(x = degree, color = label_lab)) +
  geom_freqpoly(binwidth=3, size=2) +
  facet_wrap(~Dataset, scales = "free", nrow=2) +
  scale_color_manual(values = c("Positive" = "#fe6100", "Negative" = "#785ef0")) +
  labs(x = "Degree", y = "Count", color = NULL, linetype = NULL) +
  theme_minimal(base_size = 22)+
  theme(panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.98, 0.02),
        legend.justification = c(1, 0))

ggsave("figures/fig2c_degree_density.pdf", width = 16, height = 5, dpi = 300)

#### Figure 2d: Self-interactions ####

rename_map_short <- c(
  Hippie     = "HIPPIE",
  Intact     = "IntAct",
  String     = "STRING",
  String900  = "STRING-900",
  StringDB  = "STRING-DB",
  StringExp  = "STRING-exp",
  StringText  = "STRING-text",
  PDB_dimers  = "PDB-Dimers",
  Pinder  = "PINDER"
)

self_inter <- all_interactions[, .(
  total = .N,
  self_N = sum(protein1 == protein2)
), by = c("Dataset", "label", "Split")]
self_inter[, ratio := self_N / total]
# 0-count groups can't be shown on a log scale; floor them to a half-count instead of dropping them
self_inter[, ratio_plot := ifelse(self_N == 0, 0.5 / total, ratio)]

self_inter[, label_lab := factor(ifelse(label == 1, "Positive", "Negative"), levels = c("Positive", "Negative"))]
self_inter[, Split := factor(Split, levels = c("Train", "Validation", "Test (balanced)", "Test (realistic)"))]

self_inter[Dataset %in% names(rename_map_short), Dataset := rename_map_short[Dataset]]
dataset_order <- self_inter[label_lab == "Positive", .(m = mean(ratio)), by = Dataset][order(m)]$Dataset
self_inter[, Dataset := factor(Dataset, levels = dataset_order)]

wide <- dcast(self_inter, Dataset + Split ~ label_lab, value.var = "ratio_plot")

ggplot() +
  geom_segment(data = wide, aes(x = Positive, xend = Negative, y = Dataset, yend = Dataset, group = Split),
               color = "grey85", linewidth = 0.5) +
  geom_point(data = self_inter, aes(x = ratio_plot, y = Dataset, color = label_lab, shape = Split),
             size = 3, alpha = 0.9) +
  scale_color_manual(values = c("Positive" = "#fe6100", "Negative" = "#785ef0")) +
  labs(x = "Proportion of self-interactions", y = NULL,
       color = NULL, shape = NULL,) +
  theme_minimal(base_size = 25) +
  guides(color = guide_legend(ncol = 1), shape = guide_legend(ncol=1)) +
  theme(panel.grid.minor = element_blank(), 
        legend.position = "right",
        legend.spacing.x = unit(0.1, "cm"),
        legend.margin = margin(-10, 0, 0, 0) )

ggsave("figures/fig2d_self_interactions.pdf", width = 9, height = 3.5, dpi = 300)

#### Suppl. Figure: Similarity ####

load_sim <- function(ds) {
  base <- paste0("~/PythonProjects/ppi-splitting-pipeline/results_esm2/", ds, "-Leakage-split/")
  dt <- fread(paste0(base, "similarities/all_vs_all.tsv"), sep="\t")
  colnames(dt) <- c("protein1", "protein2", "evalue", "bitscore", "pident")
  if(!ds %in% c("PDB_dimers", "Pinder")){
    dt[, protein2 := gsub(pattern="sp|", replacement="", x=protein2, fixed=TRUE)]
    dt[, protein2 := gsub(pattern="|", replacement="", x=protein2, fixed=TRUE)]
  }else if(ds == "PDB_dimers"){
    dt[, protein2 := gsub(pattern="pdb|", replacement="", x=protein2, fixed=TRUE)]
    dt[, protein2 := gsub(pattern="|", replacement="_", x=protein2, fixed=TRUE)]
  }
  dt
}

all_sim <- rbindlist(lapply(c(leakage_dirs), load_sim))
all_sim <- unique(all_sim)
all_sim_inter <- merge(all_interactions, all_sim, all.x = TRUE,
                       by.x = c("protein1", "protein2"), 
                       by.y = c("protein1", "protein2"))
all_sim_inter <- all_sim_inter[is.na(bitscore), bitscore := 0]
all_sim_inter <- all_sim_inter[is.na(evalue), evalue := 11]
all_sim_inter <- all_sim_inter[is.na(pident), pident := 0]
all_sim_inter[, label_lab := factor(ifelse(label == 1, "Positive", "Negative"), levels = c("Positive", "Negative"))]
all_sim_inter[Dataset %in% names(rename_map), Dataset := rename_map[Dataset]]
all_sim_inter[, Dataset := factor(Dataset, levels = c("PDB-Dimers",
                                                      "PINDER",
                                                      "IntAct",
                                                      "HIPPIE",
                                                      "STRING\nhuman, physical",
                                                      "STRING-\n900",
                                                      "STRING-\ndatabase",
                                                      "STRING-\nexperimental",
                                                      "STRING-\ntextmining"))]

ggplot(all_sim_inter[Split == "Train"], aes(x = pident, fill = label_lab)) +
  geom_histogram(bins=20, position='dodge') +
  facet_wrap(~Dataset, scales = "free_y", nrow=2) +
  #scale_x_log10() +
  scale_fill_manual(values = c("Positive" = "#fe6100", "Negative" = "#785ef0")) +
  labs(x = NULL, y = NULL, fill = NULL, linetype = NULL) +
  theme_minimal(base_size = 22)
ggsave("figures/fig_similarity.pdf", width = 16, height = 5, dpi = 300)

#### Figure 2e: Taxonomy ####

load_tax <- function(ds) {
  if(!ds %in% c("PDB_dimers", "Pinder")){
    base <- paste0("~/PythonProjects/ppi-splitting-pipeline/results_esm2/", ds, "-Leakage-split/")
    dt <- fread(paste0(base, "data/species.tsv"), sep="\t")
  }else if(ds == "PDB_dimers"){
    dt <- fread("~/PythonProjects/ppi-splitting-pipeline/data/pdb_dimers_species.tsv", sep="\t")
  }else{
    dt <- fread("~/PythonProjects/ppi-splitting-pipeline/data/pinder_species.tsv", sep="\t")
  }
  dt
}

all_tax <- rbindlist(lapply(c("PDB_dimers", "Pinder", "IntAct"), load_tax))
all_tax <- unique(all_tax)
all_tax_inter <- merge(all_interactions[Dataset %in% c("PDB_dimers", "Pinder", "Intact")], all_tax, by.x = "protein1", by.y = "protein_id")
setnames(all_tax_inter, "taxon_id", "taxon1")
all_tax_inter <- merge(all_tax_inter, all_tax, by.x = "protein2", by.y = "protein_id")
setnames(all_tax_inter, "taxon_id", "taxon2")
all_tax_inter <- all_tax_inter[!is.na(taxon1)]
all_tax_inter <- all_tax_inter[!is.na(taxon2)]
tax_count <- all_tax_inter[, .(
  total = .N,
  same_tax = sum(taxon1 == taxon2)
), by = c("Dataset", "label", "Split")]
tax_count[, ratio := same_tax / total]

tax_count[, label_lab := factor(ifelse(label == 1, "Positive", "Negative"), levels = c("Positive", "Negative"))]
tax_count[, Split := factor(Split, levels = c("Train", "Validation", "Test (balanced)", "Test (realistic)"))]

tax_count[Dataset %in% names(rename_map), Dataset := rename_map[Dataset]]
dataset_order <- tax_count[label_lab == "Positive", .(m = mean(ratio)), by = Dataset][order(m)]$Dataset
tax_count[, Dataset := factor(Dataset, levels = dataset_order)]

wide <- dcast(tax_count, Dataset + Split ~ label_lab, value.var = "ratio")

ggplot() +
  geom_segment(data = wide, aes(x = Positive, xend = Negative, y = Dataset, yend = Dataset, group = Split),
               color = "grey85", linewidth = 0.5) +
  geom_point(data = tax_count, aes(x = ratio, y = Dataset, color = label_lab, shape = Split),
             size = 3, alpha = 0.9) +
  scale_color_manual(values = c("Positive" = "#fe6100", "Negative" = "#785ef0")) +
  labs(x = "Proportion of same-species interactions", y = NULL,
       color = NULL, shape = NULL,) +
  theme_minimal(base_size = 16) +
  guides(color = guide_legend(ncol = 1), shape = guide_legend(ncol=2)) +
  theme(panel.grid.minor = element_blank(), 
        legend.position = "top")
ggsave("figures/fig2e_taxonomy.pdf", width = 6, height = 4, dpi = 300)

#### Figure 2f: GO terms ####

# GO-BP/MF/CC Jaccard index between the two proteins of each pair, positive vs. negative.
# Jaccard = |intersection| / |union| of GO terms (0 if the union is empty), matching
# bin/bias_analysis.py::func_relatedness. Computed via sparse protein x term indicator matrices so
# intersection/union sizes come from vectorised sparse arithmetic instead of a per-pair R loop.
#
# Pooled over Train + Validation + Test (balanced) only: those three are each ~1:1 positive:negative
# by construction, and test_balanced is a strict subset of test_realistic's pairs, so including
# test_realistic too would both double-count pairs and reintroduce its ~1:10 class imbalance.
build_go_sparse <- function(go, term_col) {
  proteins <- c(go$protein_id, "___NONE___")
  terms_split <- strsplit(go[[term_col]], ";", fixed = TRUE)
  long_protein <- rep(go$protein_id, lengths(terms_split))
  long_term <- trimws(unlist(terms_split, use.names = FALSE))
  keep <- nzchar(long_term)
  long_protein <- long_protein[keep]
  long_term <- long_term[keep]

  protein_index <- setNames(seq_along(proteins), proteins)
  term_index <- setNames(seq_along(unique(long_term)), unique(long_term))
  M <- sparseMatrix(i = protein_index[long_protein], j = term_index[long_term], x = 1,
                     dims = c(length(proteins), length(term_index)))
  list(M = M, protein_index = protein_index, sizes = Matrix::rowSums(M))
}

pooled_interactions <- all_interactions

jaccard_for_dataset <- function(ds, term_col) {
  path <- path.expand(paste0(core_path, ds, "-Leakage-split/data/go_annotations.tsv"))
  if (ds == "Pinder") {
    path <- "../../data/pinder_go_annotations.tsv"
  }else if (ds == "PDB_dimers"){
    path <- "../../data/pdbdimers_go_annotations.tsv"
  }
  go <- fread(path, select = c("protein_id", term_col))
  built <- build_go_sparse(go, term_col)

  sub <- pooled_interactions[Dataset == ds]
  idx1 <- built$protein_index[sub$protein1]
  idx1[is.na(idx1)] <- built$protein_index["___NONE___"]
  idx2 <- built$protein_index[sub$protein2]
  idx2[is.na(idx2)] <- built$protein_index["___NONE___"]

  inter <- Matrix::rowSums(built$M[idx1, , drop = FALSE] * built$M[idx2, , drop = FALSE])
  union_size <- built$sizes[idx1] + built$sizes[idx2] - inter
  jaccard <- ifelse(union_size == 0, 0, inter / union_size)
  data.table(Dataset = ds, label = sub$label, jaccard = jaccard, Split = sub$Split)
}

jac_raw <- rbindlist(lapply(list(c("go_bp", "GO-BP"), c("go_mf", "GO-MF"), c("go_cc", "GO-CC")),
                            function(x) {
                              dt <- rbindlist(lapply(leakage_dirs, jaccard_for_dataset, term_col = x[1]))
                              dt[, Category := x[2]]
                              dt
                            }))

jac_raw[, label_lab := factor(ifelse(label == 1, "Positive", "Negative"), levels = c("Positive", "Negative"))]
jac_summary <- jac_raw[, .(mean_jaccard = mean(jaccard)), by = c("Dataset", "label_lab", "Category", "Split")]

jac_summary[Dataset %in% names(rename_map_short), Dataset := rename_map_short[Dataset]]
jac_summary[, Dataset := factor(Dataset, levels = rev(c("PDB-Dimers",
                                                 "PINDER",
                                                 "HIPPIE",
                                                 "IntAct",
                                                 "STRING",
                                                 "STRING-900",
                                                 "STRING-DB",
                                                 "STRING-exp",
                                                 "STRING-text")))]
jac_summary[, Category := factor(Category, levels = c("GO-BP", "GO-MF", "GO-CC"))]


jac_wide <- dcast(jac_summary, Dataset + Category + Split ~ label_lab, value.var = "mean_jaccard")

ggplot() +
  geom_point(data = jac_summary, aes(x = mean_jaccard, y = Dataset, color = label_lab, shape = Split),
             size = 3, alpha = 0.9) +
  facet_wrap(~Category, nrow = 2) +
  scale_color_manual(values = c("Positive" = "#fe6100", "Negative" = "#785ef0")) +
  labs(x = "Jaccard index of the GO sets", y = NULL, color = NULL) +
  theme_minimal(base_size = 25) +
  theme(panel.grid.minor = element_blank(), 
        legend.position = "inside",
        legend.position.inside = c(0.9, -0.1),
        legend.justification = c(1, 0))

ggsave("figures/fig2f_go_jaccard.pdf", width = 9, height = 7, dpi = 300)



