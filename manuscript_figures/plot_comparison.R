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
library(ggforce)
library(ggtext)

core_path <- "/path/to/ppi-splitting-pipeline/results_esm2/"
path_to_results <- paste0(core_path, "multiqc_non_struc/")
path_to_struc_results <- paste0(core_path, "multiqc_structural/")
path_to_pdb_results <- paste0(core_path, "multiqc_PDB/")

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

rename_map_breaks <- c(
  Hippie     = "HIPPIE",
  Intact     = "IntAct",
  String     = "STRING\nhuman, physical",
  String900  = "STRING-900",
  StringDB  = "STRING-\ndatabase",
  StringExp  = "STRING-\nexperimental",
  StringText  = "STRING-\ntextmining",
  PDB_dimers  = "PDB-Dimers",
  Pinder  = "PINDER"
)

### Fig 3a: Dataset sizes ####

data_dirs <- c("Hippie", "Intact", "PDB_dimers", "Pinder", "String", "String900",
               "StringDB", "StringExp", "StringText")

split_files <- c(
  "Train" = "train.csv",
  "Validation" = "val.csv",
  "Test (balanced)" = "test_balanced.csv",
  "Test (realistic)" = "test_realistic.csv"
)
 
load_splits <- function(ds) {
  base <- paste0(core_path, ds)
  dt_list <- lapply(c("-Leakage-split/", "-regular/", "-ILP-split/", "-ILP-split-and-sampling/"), function(suffix){
    dt <- rbindlist(lapply(names(split_files), function(split_name) {
      d <- fread(paste0(base, suffix, split_files[[split_name]]), select = c("protein1", "protein2", "label"))
      d[, Split := split_name]
      d
    }))
    dt[, Dataset := ds]
    dt
  })
  names(dt_list) <- c("Maximally biased", "Similarity-reduced", "ILP split", "ILP split + ILP sampling")
  return_dt <- rbindlist(dt_list, idcol = "Category")
  return_dt
}

all_interactions <- rbindlist(lapply(data_dirs, load_splits))

sizes <- all_interactions[, .N, by = c("Category", "Dataset", "Split", "label")]
sizes[, Split := factor(Split, levels = c(
  "Train", "Validation", "Test (balanced)", "Test (realistic)"
))]
sizes[, Category := factor(Category, levels = c(
  "Maximally biased", "Similarity-reduced", "ILP split", "ILP split + ILP sampling"
))]
sizes[Dataset %in% names(rename_map_breaks), Dataset := rename_map_breaks[Dataset]]
sizes[, Dataset := factor(Dataset, levels = c("PDB-Dimers", 
                                              "PINDER", 
                                              "IntAct", 
                                              "HIPPIE", 
                                              "STRING\nhuman, physical", 
                                              "STRING-900",
                                              "STRING-\nexperimental",
                                              "STRING-\ndatabase",
                                              "STRING-\ntextmining"))]

sizes <- sizes[Category != "ILP split + ILP sampling"]
sizes <- sizes[Split != "Test (realistic)"]
pie_data <- sizes[label == 1] %>%
  group_by(Dataset, Category) %>%
  mutate(total = sum(N),
         end   = cumsum(N) / total * 2 * pi,
         start = lag(end, default = 0)) %>%
  ungroup() %>%
  mutate(r = sqrt(total) / max(sqrt(total)))   # area ~ N

train_labels <- pie_data %>%
  filter(Split == "Train") %>%
  distinct(Dataset, Category, N, total, r) %>%
  mutate(label = paste0(
    "n: ", scales::comma(total),
    "<br>n<sub>train</sub>: ", scales::comma(N)
  ))

ggplot(pie_data) +
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = r * 0.9,
                   start = start, end = end, fill = Split),
               color = "white", linewidth = 0.25) +
  geom_richtext(data = train_labels,
                aes(x = 0, y = 1.15, label = label),
                size = 2.6, color = "black",
                fill = NA, label.color = NA,   # no background box/border
                label.padding = unit(0, "pt"))+
  coord_fixed(clip = "off", xlim = c(-1, 1), ylim = c(-1, 1.3)) +
  facet_grid(Category ~ Dataset) +
  scale_fill_manual(values = c(
    "Train"           = "#ffb000",
    "Validation"      = "#fe6100",
    "Test (balanced)" = "#dc267f",
    "Test (realistic)"= "#785ef0"
  )) +
  labs(fill = NULL) +
  theme_void(base_size = 11) +
  theme(
    legend.position   = "bottom",
    strip.text.x      = element_text(face = "bold", size = 9, margin = margin(b = 4)),
    strip.text.y.right = element_text(face = "bold", size = 9, angle = 0, hjust = 1),
    panel.spacing     = unit(0.1, "pt"),        
    plot.margin       = margin(1, 1, 1, 1)
  )

ggsave("figures/fig3a_dataset_sizes.pdf", height=4, width=10)

#### Fig 3b: Compare similarity ####

load_sim <- function(ds) {
  base <- paste0(core_path, ds, "-Leakage-split/")
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

all_sim <- rbindlist(lapply(c(data_dirs), load_sim))
all_sim <- unique(all_sim)

prot_membership <- unique(rbind(
  all_interactions[, .(protein = protein1, Dataset, Category, Split)],
  all_interactions[, .(protein = protein2, Dataset, Category, Split)]
))

train_prots <- unique(prot_membership[Split == "Train", .(protein, Dataset, Category)])
query_prots <- prot_membership[Split != "Train"]

sim_sym <- unique(rbind(
  all_sim[, .(qry = protein1, tgt = protein2, pident, bitscore)],
  all_sim[, .(qry = protein2, tgt = protein1, pident, bitscore)]
))

hits <- merge(query_prots[, .(qry = protein, Dataset, Category, Split)],
              sim_sym, by = "qry", allow.cartesian = TRUE)

hits <- merge(hits, train_prots[, .(tgt = protein, Dataset, Category)],
              by = c("tgt", "Dataset", "Category"))

maxsim <- hits[, .(max_pident = max(pident)), by = .(Dataset, Category, Split, qry)]

maxsim <- merge(unique(query_prots[, .(qry = protein, Dataset, Category, Split)]),
                maxsim, by = c("qry", "Dataset", "Category", "Split"), all.x = TRUE)
maxsim[is.na(max_pident), max_pident := 0]

maxsim[Dataset %in% names(rename_map), Dataset := rename_map[Dataset]]
maxsim[, Dataset := factor(Dataset, levels = c("PDB-Dimers", 
                                               "PINDER", 
                                               "IntAct", 
                                               "HIPPIE", 
                                               "STRING: human, physical", 
                                               "STRING-900",
                                               "STRING-database",
                                               "STRING-experimental",
                                               "STRING-textmining"))]
maxsim[, Category := factor(Category, levels = c(
  "Maximally biased", "Similarity-reduced", "ILP split", "ILP split + ILP sampling"
))]

ggplot(maxsim[Split == "Test\n(balanced)"], aes(x = max_pident, color = Category)) +
  stat_ecdf(linewidth = 0.8) +
  geom_vline(xintercept = 40, linetype=2)+
  facet_wrap(~ Dataset, nrow=2) +
  scale_colour_manual(values = c(
    "Maximally biased"             = "#ffb000",
    "Similarity-reduced"                 = "#fe6100",
    "ILP split"              = "#dc267f",
    "ILP split + ILP sampling" = "#785ef0"
  )) +
  labs(x = "Max % identity to any training protein",
       y = "Cumulative fraction of held-out proteins",
       color = NULL) +
  theme_minimal() +
  theme(legend.position="bottom", panel.grid.minor = element_blank(),
        legend.spacing.x = unit(0.1, "cm"),
        legend.margin = margin(0, 20, 0, -10) )
ggsave("figures/fig3b_max_id.pdf", width=8.1, height=3)


### Fig 3c: Random Forest results ####
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

comparison <- rbind(test_results_balanced, test_results_real)
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
  "ILP-split-and-sampling" = "ILP split + ILP sampling"
)

comparison[Category %in% names(category_map), Category := category_map[Category]]
comparison[, Category := factor(Category, levels = rev(c(
  "Maximally biased", "Similarity-reduced", "ILP split", "ILP split + ILP sampling"
)))]

# print as LaTeX table for Appendix
comparison_dt <- comparison[order(-Dataset, -Category)]
library(kableExtra)
library(stringr)
comparison_dt[Category == "Maximally biased", Category := "Maximally\nbiased"]
comparison_dt[Category == "Similarity-reduced", Category := "Similarity-\nreduced"]
comparison_dt[Dataset == "STRING: human, physical", Dataset := "STRING\nhuman,\nphysical"]


comparison_dt %>%
  mutate(
    Dataset  = str_replace(Dataset, "-", "-\n"),
    Category = str_replace_all(Category, " \\+ ", " +\n")
  ) %>%
  mutate(
    Dataset  = linebreak(Dataset, align = "tl"),   # "tl" = top + left, not just "l"
    Category = linebreak(Category, align = "tl")
  ) %>%
  select(Dataset, Category, TestSet, AUROC, AUPRC, F1, MCC, Precision, Recall, Accuracy) %>%
  kbl(format = "latex", booktabs = TRUE, longtable = TRUE,
      digits = 2, linesep = "", escape = FALSE,
      caption = "Performance metrics across datasets, splitting strategies, and test sets.") %>%
  kable_styling(latex_options = c("repeat_header")) %>%
  collapse_rows(columns = 1:2, latex_hline = "major", valign = "top") %>%   # valign = "top" here too
  column_spec(1, width = "2.0cm") %>%
  column_spec(2, width = "2.2cm") %>%
  save_kable(file = "figures/rf_table.tex")

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
    "Maximally biased" = "#ffb000",
    "Similarity-reduced" = "#fe6100",
    "ILP split" = "#dc267f",
    "ILP split + ILP sampling" = "#785ef0"
  )) +
  scale_x_reverse() +                                    
  coord_cartesian(xlim = c(1.1, 0.48), expand = FALSE)+
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.y  = element_blank(),   
        axis.ticks.y = element_blank(),   
        axis.title.y = element_blank(),   
        plot.margin = margin(5.5, 2, 5.5, 5.5)) +
  ggtitle("Balanced")

p_realistic <- ggplot(comparison %>% filter(TestSet == "realistic"), 
                      aes(y = Dataset, x = AUPRC, fill = Category)) +
  geom_col(position = dodge_width) +
  geom_text(aes(label = sprintf("%.2f", AUPRC)),
            position = dodge_width,
            hjust = -0.1,
            size = 3.2) +
  geom_vline(xintercept = 1/11, linetype = 2) +
  scale_fill_manual(values = c(
    "Maximally biased" = "#ffb000",
    "Similarity-reduced" = "#fe6100",
    "ILP split" = "#dc267f",
    "ILP split + ILP sampling" = "#785ef0"
  )) +
  coord_cartesian(xlim = c(0.08, 1.1), expand = FALSE) +
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
ggsave("figures/fig3c_rf_comparison.pdf", height=5, width=8)



