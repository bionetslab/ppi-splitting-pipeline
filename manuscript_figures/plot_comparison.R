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

path_to_robustness <- "/path/to/ppi-splitting-pipeline/results_robustness/"
path_to_hcn <- "/path/to/ppi-splitting-pipeline/results_robustness_hcn/"

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

category_map <- c(
  "Leakage-split" = "Maximally biased",
  "regular" = "Similarity-reduced",
  "ILP-split" = "ILP split",
  "ILP-split-and-sampling" = "ILP split + ILP sampling",
  "ILP-hcn" = "ILP split + ILP sampling (HRNI)"
)

pop_sd <- function(x) sqrt(mean((x - mean(x))^2))  # population SD, as in bin/robustness_summary.py

split_files <- c(
  "Train" = "train.csv",
  "Validation" = "val.csv",
  "Test (balanced)" = "test_balanced.csv"
)

# Positives only depend on the split seed, so one sampling seed (42) per split seed suffices
load_splits <- function(run_dir) {
  rbindlist(lapply(Sys.glob(paste0(run_dir, "/split*")), function(split_dir) {
    sampling_dir <- Sys.glob(paste0(split_dir, "/sampling42"))
    rbindlist(lapply(names(split_files), function(split_name) {
      d <- fread(paste0(sampling_dir, "/", split_files[[split_name]]), select = "label")
      data.table(Split = split_name, N = sum(d$label == 1))
    }))[, split_seed := basename(split_dir)]
  }))[, ID := basename(dirname(run_dir))]
}

sizes_per_seed <- rbindlist(lapply(Sys.glob(paste0(path_to_robustness, "*/robustness")), load_splits))
sizes_per_seed[, `:=`(Dataset = sub("-.*", "", ID), Category = category_map[sub("^[^-]*-", "", ID)])]
sizes_per_seed[, total := sum(N), by = .(ID, split_seed)]

# Pies show the mean over split seeds
sizes <- sizes_per_seed[, .(N = mean(N), N_sd = pop_sd(N), total_sd = pop_sd(total)),
                        by = .(Category, Dataset, Split)]
sizes[, Split := factor(Split, levels = names(split_files))]
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
pie_data <- sizes %>%
  group_by(Dataset, Category) %>%
  mutate(total = sum(N),
         end   = cumsum(N) / total * 2 * pi,
         start = lag(end, default = 0)) %>%
  ungroup() %>%
  mutate(r = sqrt(total) / max(sqrt(total)))   # area ~ N

short_num <- scales::label_number(accuracy = 0.1, scale_cut = scales::cut_short_scale(), drop0trailing = TRUE)
train_labels <- pie_data %>%
  filter(Split == "Train") %>%
  distinct(Dataset, Category, N, N_sd, total, total_sd, r) %>%
  mutate(label = paste0(
    "n: ", short_num(total), " ± ", short_num(total_sd),
    "<br>n<sub>train</sub>: ", short_num(N), " ± ", short_num(N_sd)
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

# Held-out and training proteins per split seed; the protein sets do not depend on the sampling seed
load_proteins <- function(run_dir) {
  rbindlist(lapply(Sys.glob(paste0(run_dir, "/split*")), function(split_dir) {
    sampling_dir <- Sys.glob(paste0(split_dir, "/sampling*"))[1]
    rbindlist(lapply(c(Train = "train.csv", Test = "test_balanced.csv"), function(f) {
      d <- fread(paste0(sampling_dir, "/", f), select = c("protein1", "protein2"))
      data.table(protein = unique(c(d$protein1, d$protein2)))
    }), idcol = "Split")[, split_seed := basename(split_dir)]
  }))[, ID := basename(dirname(run_dir))]
}

prot_membership <- rbindlist(lapply(Sys.glob(paste0(path_to_robustness, "*/robustness")), load_proteins))
prot_membership[, `:=`(Dataset = sub("-.*", "", ID), Category = category_map[sub("^[^-]*-", "", ID)])]

load_sim <- function(ds) {
  # all_vs_all.tsv is identical for every splitting strategy of a dataset
  dt <- fread(Sys.glob(paste0(path_to_robustness, ds, "-*/similarities/all_vs_all.tsv"))[1], sep="\t")
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

# Max % identity of each test protein to any training protein, per splitting strategy and split seed
max_identity <- function(ds) {
  sim <- load_sim(ds)
  sim <- unique(rbind(sim[, .(qry = protein1, tgt = protein2, pident)],
                      sim[, .(qry = protein2, tgt = protein1, pident)]))
  prots <- prot_membership[Dataset == ds]
  prots[Split == "Test", {
    train <- prots[Split == "Train" & Category == .BY$Category & split_seed == .BY$split_seed, protein]
    res <- sim[qry %in% protein & tgt %in% train, .(max_pident = max(pident)), by = qry]
    res <- res[data.table(qry = protein), on = "qry"]
    res[is.na(max_pident), max_pident := 0]
    res
  }, by = .(Category, split_seed)][, Dataset := ds]
}

maxsim <- rbindlist(lapply(unique(prot_membership$Dataset), max_identity))

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

# ILP split + ILP sampling lies on top of ILP split
ggplot(maxsim[Category != "ILP split + ILP sampling"],
       aes(x = max_pident, color = Category, group = interaction(Category, split_seed))) +
  stat_ecdf(linewidth = 0.4) +  # one line per split seed
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
dataset_levels <- rev(c("PDB-Dimers", 
                        "PINDER", 
                        "IntAct", 
                        "HIPPIE", 
                        "STRING: human, physical", 
                        "STRING-900",
                        "STRING-experimental",
                        "STRING-database",
                        "STRING-textmining"))

metrics <- c("AUROC", "AUPRC", "F1", "MCC", "Precision", "Recall", "Accuracy")

# Seed robustness runs (split x sampling x training seeds); datasets appear once their runs are in path_to_robustness
robustness_runs <- rbindlist(lapply(Sys.glob(paste0(c(path_to_robustness, path_to_hcn), "*/robustness/robustness_runs.tsv")),
                                    fread, sep = "\t"))
robustness_runs <- robustness_runs[ID != "Intact-ILP-hcn"]  # IntAct is multi-species, the HRNI network is human-only
robustness_runs[, `:=`(
  Dataset  = factor(rename_map[sub("-.*", "", ID)], levels = dataset_levels),
  Category = factor(category_map[sub("^[^-]*-", "", ID)], levels = rev(unname(category_map))),
  TestSet  = sub("test_", "", test_set)
)]

# print as LaTeX table for Appendix: mean ± SD over seed runs
comparison_dt <- robustness_runs[, lapply(.SD, function(x) sprintf("%.2f $\\pm$ %.2f", mean(x), pop_sd(x))),
                                 by = .(Dataset, Category, TestSet), .SDcols = metrics][order(-Dataset, -Category, TestSet)]
library(kableExtra)
library(stringr)
comparison_dt[Category == "Maximally biased", Category := "Maximally\nbiased"]
comparison_dt[Category == "Similarity-reduced", Category := "Similarity-\nreduced"]
comparison_dt[Category == "ILP split + ILP sampling (HRNI)", Category := "ILP split + sampling (HRNI)"]  # shorter, fits the column

# Dataset becomes a group-header row (\multicolumn{9}{l}{Dataset: \textbf{...}}) instead of a column
dataset_runs <- rle(as.character(comparison_dt$Dataset))
dataset_index <- setNames(dataset_runs$lengths, paste0("Dataset: \\textbf{", dataset_runs$values, "}"))

comparison_dt %>%
  mutate(Category = linebreak(str_replace_all(Category, " \\+ ", " +\n"), align = "tl")) %>%  # "tl" = top + left
  select(Category, TestSet, AUROC, AUPRC, F1, MCC, Precision, Recall, Accuracy) %>%
  kbl(format = "latex", booktabs = TRUE, longtable = TRUE,
      digits = 2, linesep = "", escape = FALSE,
      col.names = c("Category", "Test set", "AUROC", "AUPRC", "F1", "MCC", "Precision", "Recall", "Accuracy"),
      caption = "Performance metrics across datasets, splitting strategies, and test sets.") %>%
  kable_styling(latex_options = c("repeat_header")) %>%
  collapse_rows(columns = 1, latex_hline = "none", valign = "top") %>%  # must run before pack_rows, else it mangles \textbf
  column_spec(1, width = "1.8cm") %>%  # after collapse_rows: the other order drops a } from \makecell cells
  pack_rows(index = dataset_index, escape = FALSE, bold = FALSE, indent = FALSE, latex_gap_space = "0em") %>%
  # \midrule above and below each Dataset row; pack_rows only offers \hline, which clashes with booktabs
  gsub("\\\\addlinespace\\[0em\\]\n", "\\\\midrule\n", .) %>%
  gsub("(\\\\multicolumn\\{[0-9]+\\}\\{l\\}\\{Dataset:[^\n]*\n)", "\\1\\\\midrule\n", .) %>%
  sub("\\\\endlastfoot\n\\\\midrule\n", "\\\\endlastfoot\n", .) %>%  # first group already sits under the header rule
  paste0("{\\small\\setlength{\\tabcolsep}{3pt}\n", ., "}") %>%  # braces keep font/spacing local to this table
  writeLines("figures/rf_table.tex")

# Bars show the mean over seed runs
plot_data <- robustness_runs[, .(
  AUPRC      = mean(AUPRC),
  AUPRC_sd   = pop_sd(AUPRC),
  AUPRC_max  = max(AUPRC)
), by = .(Dataset, Category, TestSet)]
plot_data <- plot_data[Category != "ILP split + ILP sampling (HRNI)"]  # HRNI goes in the table only
plot_data[, label_x := pmax(AUPRC + AUPRC_sd, AUPRC_max)]  # labels beyond error bars and dots

dodge_width <- position_dodge(width = 0.9)

robustness_layers <- function(panel) list(
  # group = Category keeps the dots dodged like the bars instead of splitting them per split seed
  #geom_point(data = robustness_runs %>% filter(TestSet == panel),
  #           aes(shape = factor(split_seed), group = Category),
  #           position = dodge_width, size = 1.3, alpha = 0.5, colour = "grey20",
  #           show.legend = c(shape = TRUE, fill = FALSE)),
  geom_errorbar(aes(xmin = AUPRC - AUPRC_sd, xmax = AUPRC + AUPRC_sd),
                position = dodge_width, width = 0.5, linewidth = 0.4, show.legend = FALSE)
)

p_balanced <- ggplot(plot_data %>% filter(TestSet == "balanced"),
                     aes(y = Dataset, x = AUPRC, fill = Category)) +
  geom_col(position = dodge_width) +
  robustness_layers("balanced") +
  geom_text(aes(x = label_x, label = sprintf("%.2f ± %.2f", AUPRC, AUPRC_sd)),
            position = dodge_width,
            hjust = 1.05,
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

p_realistic <- ggplot(plot_data %>% filter(TestSet == "realistic"), 
                      aes(y = Dataset, x = AUPRC, fill = Category)) +
  geom_col(position = dodge_width) +
  robustness_layers("realistic") +
  geom_text(aes(x = label_x, label = sprintf("%.2f ± %.2f", AUPRC, AUPRC_sd)),
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
legend_source <- p_realistic + labs(fill = NULL, shape = "Split seed")+
  guides(shape = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  theme(legend.position="bottom")
shared_legend <- cowplot::get_legend(legend_source)

# Combine the two panels side by side
plots_row <- plot_grid(p_balanced, p_realistic, nrow = 1,
                       align = "h", axis = "tb",
                       rel_widths = c(1, 1))

# Stack the row of plots on top of the shared legend
final_plot <- plot_grid(plots_row, shared_legend, ncol = 1, rel_heights = c(1, 0.1))
final_plot
ggsave("figures/fig3c_rf_comparison.pdf", height=5.5, width=9)



