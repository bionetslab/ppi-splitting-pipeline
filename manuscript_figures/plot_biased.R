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

path_robustness_bias <- "/path/to/project/ppi-splitting-pipeline/results_robustness_bias/"
path_to_robustness <- "/path/to/project/ppi-splitting-pipeline/results_robustness/"

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

all_biases <- lapply(list.files(path_robustness_bias, pattern="*mqc.tsv", full.names=TRUE, recursive=TRUE), fread)
names(all_biases) <- list.files(path_robustness_bias, pattern="*mqc.tsv", full.names=TRUE, recursive=TRUE)

all_biases <- rbindlist(all_biases, idcol="filename")
all_biases[, filename := gsub(paste0(path_robustness_bias, "/"), "", filename, fixed=TRUE)]
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

metrics <- c("AUROC", "AUPRC", "F1", "MCC", "Precision", "Recall", "Accuracy")

# Means over all seed combinations of the robustness runs
robustness_runs <- rbindlist(lapply(Sys.glob(paste0(path_to_robustness, "*-Leakage-split/robustness/robustness_runs.tsv")),
                                    fread, sep = "\t"))
only_leakage <- robustness_runs[, lapply(.SD, mean),
                                by = .(Dataset = sub("-.*", "", ID), TestSet = sub("test_", "", test_set)),
                                .SDcols = metrics]
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

only_leakage <- only_leakage[order(Dataset, TestSet)]

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

# Training degrees for every (split, sampling) seed combination of the robustness runs; the training seed plays no role
load_train_degrees <- function(f) {
  d <- fread(f, select = c("protein1", "protein2", "label"))
  # a self-interaction (protein1 == protein2) counts once
  long <- rbind(d[, .(protein = protein1, label)], d[protein1 != protein2, .(protein = protein2, label)])
  long[, .(degree = .N), by = .(label, protein)][, `:=`(
    Dataset = sub("-.*", "", basename(dirname(dirname(dirname(dirname(f)))))),
    seeds   = paste(basename(dirname(dirname(f))), basename(dirname(f)))
  )]
}

degree_dt <- rbindlist(lapply(Sys.glob(paste0(path_to_robustness, "*-Leakage-split/robustness/split*/sampling*/train.csv")),
                              load_train_degrees))

degree_dt[, label_lab := factor(ifelse(label == 1, "Positive", "Negative"), levels = c("Positive", "Negative"))]

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

ggplot(degree_dt[degree <=100], aes(x = degree, color = label_lab, group = interaction(label_lab, seeds))) +
  geom_freqpoly(binwidth=3, linewidth=0.5) +  # one line per (split, sampling) seed combination
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

# Every split file of every (split, sampling) seed combination of the Maximally biased robustness runs
split_names <- c(train = "Train", val = "Validation", test_balanced = "Test (balanced)", test_realistic = "Test (realistic)")
leakage_runs <- data.table(file = Sys.glob(paste0(path_to_robustness, "*-Leakage-split/robustness/split*/sampling*/*.csv")))
leakage_runs[, `:=`(
  Dataset = sub("-.*", "", basename(dirname(dirname(dirname(dirname(file)))))),
  seeds   = paste(basename(dirname(dirname(file))), basename(dirname(file))),
  Split   = split_names[sub(".csv", "", basename(file), fixed = TRUE)]
)]

# Applies stat(pairs, dataset) to each run file, one file at a time to keep memory bounded
per_run <- function(runs, stat) {
  rbindlist(lapply(seq_len(nrow(runs)), function(i) {
    d <- fread(runs$file[i], select = c("protein1", "protein2", "label"))
    stat(d, runs$Dataset[i])[, `:=`(Dataset = runs$Dataset[i], Split = runs$Split[i], seeds = runs$seeds[i])]
  }))
}

self_runs <- per_run(leakage_runs, function(d, ds) d[, .(total = .N, ratio = mean(protein1 == protein2)), by = label])
# One row (= one point) per (split, sampling) seed combination
self_inter <- self_runs

# Offsets the splits vertically within each dataset row
split_dodge <- position_dodge(width = 0.7)

self_inter[, label_lab := factor(ifelse(label == 1, "Positive", "Negative"), levels = c("Positive", "Negative"))]
self_inter[, Split := factor(Split, levels = c("Train", "Validation", "Test (balanced)", "Test (realistic)"))]

self_inter[Dataset %in% names(rename_map_short), Dataset := rename_map_short[Dataset]]
dataset_order <- self_inter[label_lab == "Positive", .(m = mean(ratio)), by = Dataset][order(m)]$Dataset
self_inter[, Dataset := factor(Dataset, levels = dataset_order)]


ggplot() +
  geom_point(data = self_inter, aes(x = ratio, y = Dataset, color = label_lab, shape = Split, group = Split),
             size = 2) +
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
  dt <- fread(paste0(path_to_robustness, ds, "-Leakage-split/similarities/all_vs_all.tsv"), sep="\t")
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

# % identity of the training pairs of one seed combination (split 1234, sampling 42; training seed irrelevant here);
# pairs without a BLAST hit count as 0 %
sim_pairs <- rbindlist(lapply(unique(leakage_runs$Dataset), function(ds) {
  sim <- unique(load_sim(ds)[, .(protein1, protein2, pident)])
  per_run(leakage_runs[Dataset == ds & Split == "Train" & seeds == "split1234 sampling42"], function(d, ds) {
    d <- merge(d, sim, by = c("protein1", "protein2"), all.x = TRUE)
    d[is.na(pident), pident := 0]
    d[, .(label, pident)]
  })
}))
sim_pairs[, label_lab := factor(ifelse(label == 1, "Positive", "Negative"), levels = c("Positive", "Negative"))]
sim_pairs[Dataset %in% names(rename_map), Dataset := rename_map[Dataset]]
sim_pairs[, Dataset := factor(Dataset, levels = c("PDB-Dimers",
                                                      "PINDER",
                                                      "IntAct",
                                                      "HIPPIE",
                                                      "STRING\nhuman, physical",
                                                      "STRING-\n900",
                                                      "STRING-\ndatabase",
                                                      "STRING-\nexperimental",
                                                      "STRING-\ntextmining"))]

ggplot(sim_pairs, aes(x = pident, fill = label_lab)) +
  geom_histogram(bins=20, position='dodge') +
  scale_y_continuous(trans = "log1p", breaks = c(0, 100, 1000, 10000, 100000))+
  facet_wrap(~Dataset, scales = "free_y", nrow=2) +
  scale_fill_manual(values = c("Positive" = "#fe6100", "Negative" = "#785ef0")) +
  labs(x = NULL, y = NULL, fill = NULL, linetype = NULL) +
  theme_minimal(base_size = 22)
ggsave("figures/FigureA1.pdf", width = 16, height = 5, dpi = 300)

#### Figure 2e: Taxonomy ####

load_tax <- function(ds) {
  if(!ds %in% c("PDB_dimers", "Pinder")){
    dt <- fread(paste0(path_to_robustness, ds, "-Leakage-split/data/species.tsv"), sep="\t")
  }else if(ds == "PDB_dimers"){
    dt <- fread("/path/to/project/ppi-splitting-pipeline/data/data/pdb_dimers_species.tsv", sep="\t")
  }else{
    dt <- fread("/path/to/project/ppi-splitting-pipeline/data/data/pinder_species.tsv", sep="\t")
  }
  dt
}

tax_datasets <- intersect(c("PDB_dimers", "Pinder", "Intact"), leakage_runs$Dataset)
all_tax <- rbindlist(lapply(tax_datasets, load_tax))
all_tax <- unique(all_tax)

tax_runs <- per_run(leakage_runs[Dataset %in% tax_datasets], function(d, ds) {
  d <- merge(d, all_tax, by.x = "protein1", by.y = "protein_id")
  setnames(d, "taxon_id", "taxon1")
  d <- merge(d, all_tax, by.x = "protein2", by.y = "protein_id")
  setnames(d, "taxon_id", "taxon2")
  d[!is.na(taxon1) & !is.na(taxon2), .(ratio = mean(taxon1 == taxon2)), by = label]
})
# One row (= one point) per (split, sampling) seed combination
tax_count <- tax_runs

tax_count[, label_lab := factor(ifelse(label == 1, "Positive", "Negative"), levels = c("Positive", "Negative"))]
tax_count[, Split := factor(Split, levels = c("Train", "Validation", "Test (balanced)", "Test (realistic)"))]

tax_count[Dataset %in% names(rename_map), Dataset := rename_map[Dataset]]
dataset_order <- tax_count[label_lab == "Positive", .(m = mean(ratio)), by = Dataset][order(m)]$Dataset
tax_count[, Dataset := factor(Dataset, levels = dataset_order)]


ggplot() +
  geom_point(data = tax_count, aes(x = ratio, y = Dataset, color = label_lab, shape = Split, group = Split),
             size = 2) +
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

go_categories <- c("GO-BP" = "go_bp", "GO-MF" = "go_mf", "GO-CC" = "go_cc")

jaccard_for_dataset <- function(ds) {
  path <- path.expand(paste0(path_to_robustness, ds, "-Leakage-split/data/go_annotations.tsv"))
  if (ds == "Pinder") {
    path <- "/path/to/project/ppi-splitting-pipeline/data/data/pinder_go_annotations.tsv"
  }else if (ds == "PDB_dimers"){
    path <- "/path/to/project/ppi-splitting-pipeline/data/data/pdb_dimers_go_annotations.tsv"
  }
  go <- fread(path, select = c("protein_id", unname(go_categories)))
  built <- lapply(go_categories, function(term_col) build_go_sparse(go, term_col))

  per_run(leakage_runs[Dataset == ds], function(d, ds) {
    rbindlist(lapply(names(built), function(category) {
      b <- built[[category]]
      idx1 <- b$protein_index[d$protein1]
      idx1[is.na(idx1)] <- b$protein_index["___NONE___"]
      idx2 <- b$protein_index[d$protein2]
      idx2[is.na(idx2)] <- b$protein_index["___NONE___"]

      inter <- Matrix::rowSums(b$M[idx1, , drop = FALSE] * b$M[idx2, , drop = FALSE])
      union_size <- b$sizes[idx1] + b$sizes[idx2] - inter
      jaccard <- ifelse(union_size == 0, 0, inter / union_size)
      data.table(label = d$label, jaccard = jaccard)[, .(mean_jaccard = mean(jaccard)), by = label][, Category := category]
    }))
  })
}

jac_runs <- rbindlist(lapply(unique(leakage_runs$Dataset), jaccard_for_dataset))

jac_runs[, label_lab := factor(ifelse(label == 1, "Positive", "Negative"), levels = c("Positive", "Negative"))]
# One row (= one point) per (split, sampling) seed combination
jac_summary <- jac_runs

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
jac_summary[, Split := factor(Split, levels = split_names)]


jac_wide <- dcast(jac_summary, Dataset + Category + Split ~ label_lab, value.var = "mean_jaccard", fun.aggregate = mean)  # mean over seed combinations

ggplot() +
  geom_point(data = jac_summary, aes(x = mean_jaccard, y = Dataset, color = label_lab, shape = Split, group = Split),
             size = 2) +
  facet_wrap(~Category, nrow = 2) +
  scale_color_manual(values = c("Positive" = "#fe6100", "Negative" = "#785ef0")) +
  labs(x = "Jaccard index of the GO sets", y = NULL, color = NULL) +
  theme_minimal(base_size = 25) +
  theme(panel.grid.minor = element_blank(), 
        legend.position = "inside",
        legend.position.inside = c(0.9, -0.1),
        legend.justification = c(1, 0))

ggsave("figures/fig2f_go_jaccard.pdf", width = 9, height = 7, dpi = 300)



