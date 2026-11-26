library(data.table)
library(ggplot2)
library(stringr)

path_to_robustness <- "/path/to/ppi-splitting-pipeline/results_robustness/"
path_to_input <- "/path/to/ppi-splitting-pipeline/data/data/"

pop_sd <- function(x) sqrt(mean((x - mean(x))^2))  # population SD, as in bin/robustness_summary.py

#### Plot data loss from CD-Hit ####

input_ppis <- c(Hippie = "HIPPIE-current.csv", Intact = "intact.csv", String = "string.csv", String900 = "string_900.csv",
                StringDB = "string_database.csv", StringExp = "string_experimental.csv", StringText = "string_textmining.csv",
                PDB_dimers = "pdb_interactions.csv", Pinder = "pinder.csv")

# Rebuilds REMOVE_REDUNDANT's split bar counts for one (run, split seed); robustness runs don't publish them.
# Each KaHIP cluster's split is read off the final splits (majority of its surviving proteins), giving the
# pre-CD-HIT split sizes: KaHIP/ILP discards = input - pre-CD-HIT, CD-HIT discards = pre-CD-HIT - kept.
# Positives only depend on the split seed, so sampling seed 42 suffices.
split_bar_counts <- function(split_dir, ppis, node_map) {
  part <- data.table(protein_id = node_map$protein_id,
                     cluster = scan(paste0(split_dir, "/partitioned_proteome.txt"), quiet = TRUE)[node_map$node_id])
  kept <- rbindlist(lapply(c(Train = "train", Val = "val", Test = "test_balanced"), function(f) {
    d <- fread(paste0(split_dir, "/sampling42/", f, ".csv"))[label == 1]
    data.table(n = nrow(d), protein_id = unique(c(d$protein1, d$protein2)))
  }), idcol = "Sample")
  cluster_split <- merge(kept, part, by = "protein_id")[, .N, by = .(cluster, Sample)][order(-N)][, .SD[1], by = cluster]
  part <- merge(part, cluster_split[, .(cluster, Sample)], by = "cluster")
  s1 <- part$Sample[match(ppis$protein1, part$protein_id)]
  s2 <- part$Sample[match(ppis$protein2, part$protein_id)]
  n_pre <- sum(!is.na(s1) & s1 == s2)
  n_kept <- unique(kept[, .(Sample, n)])
  data.table(Sample = c(n_kept$Sample, "Discarded"),
             Kept = c(n_kept$n, 0),
             `Discarded (KaHIP/ILP)` = c(0, 0, 0, nrow(ppis) - n_pre),
             `Discarded (CD-HIT-2D)` = c(0, 0, 0, n_pre - sum(n_kept$n)))
}

cdhit_runs <- rbindlist(lapply(Sys.glob(paste0(path_to_robustness, "*-", c("regular", "ILP-split"))), function(run_dir) {
  id <- basename(run_dir)
  ppis <- fread(paste0(path_to_input, input_ppis[[sub("-.*", "", id)]]), select = c("protein1", "protein2"))
  # unique unordered pairs, as bin/utils.py::read_ppis
  ppis <- unique(ppis[, .(protein1 = trimws(pmin(protein1, protein2)), protein2 = trimws(pmax(protein1, protein2)))])
  node_map <- fread(paste0(run_dir, "/similarities/node_mapping.tsv"))
  rbindlist(lapply(Sys.glob(paste0(run_dir, "/robustness/split*")), function(split_dir)
    split_bar_counts(split_dir, ppis, node_map)[, split_seed := basename(split_dir)]))[, ID := id]
}))
cdhit_runs[, `:=`(Dataset = sub("-.*", "", ID), Category = sub("^[^-]*-", "", ID))]

# Bars show the mean over split seeds, error bars the SD
cdhit_plt <- melt(cdhit_runs, id.vars = c("ID", "split_seed", "Sample", "Dataset", "Category"))
cdhit_plt <- cdhit_plt[, .(value = mean(value), sd = pop_sd(value)), by = .(Sample, Dataset, Category, variable)][value > 0]

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

cdhit_plt[Dataset %in% names(rename_map_breaks), Dataset := rename_map_breaks[Dataset]]
cdhit_plt[, Dataset := factor(Dataset, levels = c("PDB-Dimers",
                                                  "PINDER",
                                                  "IntAct",
                                                  "HIPPIE",
                                                  "STRING\nhuman, physical",
                                                  "STRING-900",
                                                  "STRING-\nexperimental",
                                                  "STRING-\ndatabase",
                                                  "STRING-\ntextmining"))]

cdhit_plt[, Sample := factor(Sample, levels = c("Train", "Val", "Test", "Discarded"))]
cdhit_plt[, Category := factor(c(regular = "Similarity-reduced", `ILP-split` = "ILP split")[Category],
                               levels = c("Similarity-reduced", "ILP split"))]

# Error bars sit at the end of each stacked segment; position_stack puts the first fill level outermost
setorder(cdhit_plt, Dataset, Category, Sample, -variable)
cdhit_plt[, end := cumsum(value), by = .(Dataset, Category, Sample)]

ggplot(cdhit_plt, aes(x = value, y = Sample, fill = variable))+
  geom_col()+
  geom_errorbar(aes(xmin = end - sd, xmax = end + sd), width = 0.3, linewidth = 0.4)+
  facet_wrap(Dataset~Category, scales = "free", ncol=4)+
  scale_fill_manual(values = c(
    "Kept"             = "#fe6100",
    "Discarded (KaHIP/ILP)"                 = "#dc267f",
    "Discarded (CD-HIT-2D)"              = "#ffb000"
  )) +
  scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                       breaks = scales::pretty_breaks(n = 2)) +
  theme_minimal(base_size=15)+
  theme(strip.text.x = element_text(face = "bold", size = 15, margin = margin(b = 4)),
        legend.position = "bottom"
        )+
  labs(fill = NULL, x = NULL, y = NULL)
ggsave("figures/FigureA2.pdf", width=12, height=15)
