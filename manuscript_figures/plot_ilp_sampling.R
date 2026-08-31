library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

core_path <- "/path/to/ppi-splitting-pipeline/results_server/"

### Fig 4a: Bias analysis results ####

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
all_biases <- all_biases[Category != "ILP-hcn"]
all_biases <- all_biases[Bias != "topology_shortcut"]

all_biases[Dataset %in% names(rename_map), Dataset := rename_map[Dataset]]
all_biases[, Dataset := factor(Dataset, levels = c("PDB-Dimers", 
                                                   "PINDER", 
                                                   "IntAct", 
                                                   "HIPPIE", 
                                                   "STRING\nhuman, physical", 
                                                   "STRING-\n900",
                                                   "STRING-\nexperimental",
                                                   "STRING-\ndatabase",
                                                   "STRING-\ntextmining"))]

rename_split_map <- c(
  train = "Train",
  val = "Validation",
  test_balanced = "Test (balanced)",
  test_realistic = "Test (realistic)"
)
all_biases[Split %in% names(rename_split_map), Split := rename_split_map[Split]]
all_biases[, Split := factor(Split, levels = c("Train", "Validation", "Test (balanced)", "Test (realistic)"))]

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
all_biases[Bias %in% names(rename_bias_map), Bias := rename_bias_map[Bias]]
all_biases[, Bias := factor(Bias, levels = c("Self interactions", 
                                             "Same Species", 
                                             "Topology",
                                             "Sequence Similarity",
                                             "Embedding Similarity",
                                             "GO-BP",
                                             "GO-MF",
                                             "GO-CC"))]

category_map <- c(
  "Leakage-split" = "Maximally biased",
  "regular" = "Similarity-reduced",
  "ILP-split" = "ILP-based split",
  "ILP-split-and-sampling" = "ILP-based split and sampling"
)

all_biases[Category %in% names(category_map), Category := category_map[Category]]
all_biases[, Category := factor(Category, levels = c(
  "Maximally biased", "Similarity-reduced", "ILP-based split", "ILP-based split and sampling"
))]

all_biases[`Detectability (Spearman ρ)` == "NaN", `Detectability (Spearman ρ)` := NA]

ggplot(all_biases[Category %in% c("ILP-based split and sampling", "ILP-based split")], 
       aes(x = Category, y = `NMI(A;Y)`, group = Bias, color = Bias)) +
  geom_line(linewidth = 0.8) +
  scale_x_discrete(labels = c("ILP-based split" = "ILP Split",
                              "ILP-based split and sampling" = "ILP Split +sampling"))+
  geom_point(size = 2.5) +
  scale_color_manual(values=c("Self interactions" = "#ffb000", 
                              "Sequence Similarity"="#fe6100", 
                              "Same Species"="#fc92d7",
                              "Embedding Similarity" = "#c701ff",
                              "GO-BP" = "#648fff",
                              "GO-MF" = "#003a7d",
                              "GO-CC" = "#bdd373"))+
  facet_grid(Split ~ Dataset) +
  guides(color = guide_legend(nrow = 1)) +
  labs(y = "NMI(a; y)", x = NULL, color = NULL) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "top",
        strip.text.x      = element_text(face = "bold", size = 15, margin = margin(b = 4)),
        strip.text.y      = element_text(face = "bold", size = 15, margin = margin(b = 4)),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1))

ggsave("figures/fig4a_bias_analysis_comparison.pdf", width = 15, height = 8, dpi = 300)

### Fig 4b: Bias objective terms results ###

path_to_results <- paste0(core_path, "multiqc_non_struc/")
path_to_struc_results <- paste0(core_path, "multiqc_structural/")
path_to_pdb_results <- paste0(core_path, "multiqc_PDB/")

residuals <- fread(paste0(path_to_results, "multiqc_report_data/multiqc_neg_ilp_diagnostics_table.txt"), sep="\t")
residuals <- residuals[, c("ID", "split", "obj_value", "bias_deg_term", "bias_tax_term", "bias_self_term", "bias_jac_term", "status")]
res_struc <- fread(paste0(path_to_struc_results, "multiqc_report_data/multiqc_neg_ilp_diagnostics_table.txt"), sep="\t")
res_struc <- res_struc[!startsWith(Sample, "PDB")]
res_struc <- res_struc[, c("ID", "split", "obj_value", "bias_deg_term", "bias_tax_term", "bias_self_term", "bias_jac_term", "status")]
res_pdb <- fread(paste0(path_to_pdb_results, "multiqc_report_data/multiqc_neg_ilp_diagnostics_table.txt"), sep="\t")
res_pdb <- res_pdb[, c("ID", "split", "obj_value", "bias_deg_term", "bias_tax_term", "bias_self_term", "bias_jac_term", "status")]
residuals <- rbind(residuals, res_struc, res_pdb)
residuals$Dataset <- sapply(stringr::str_split(residuals$ID, "-", n = 2),`[`, 1)

res_long <- melt(residuals, id.vars = c("Dataset", "ID", "split", "obj_value", "status"))
res_long[split %in% names(rename_split_map), split := rename_split_map[split]]
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
  "Taxonomy",
  "GO-BP Jaccard index"
))]

res_long[split == "Test (balanced)", split := "Test (bal.)"]
res_long[split == "Validation", split := "Val."]
res_long[, split := factor(split, levels = c("Train", "Val.", "Test (bal.)", "Test (realistic)"))]

ggplot(res_long, aes(x = split, y = value, fill = variable))+
  geom_col()+
  facet_wrap(~Dataset, nrow=2)+
  scale_fill_manual(values = c(
    "Degree"              = "#ffb000",
    "Self interactions"             = "#fe6100",
    "Taxonomy"                 = "#dc267f",
    "GO-BP Jaccard index" = "#785ef0"
  )) +
  theme_minimal(base_size=15)+
  labs(fill = NULL, y = "Objective value", x = NULL)+
  theme(legend.position = "inside",
        legend.position.inside = c(0.98, 0.02),
        legend.justification = c(1, 0),
        panel.grid.minor = element_blank(),
        strip.text.x      = element_text(face = "bold", size = 15, margin = margin(b = 4)),
        )
ggsave("figures/fig4b_residuals.pdf", width=12, height=5, dpi=300)


