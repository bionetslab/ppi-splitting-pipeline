library(data.table)
library(ggplot2)
library(stringr)

core_path <- "/path/to/ppi-splitting-pipeline/results_server/"
path_to_results <- paste0(core_path, "multiqc_non_struc/")
path_to_struc_results <- paste0(core_path, "multiqc_structural/")
path_to_pdb_results <- paste0(core_path, "multiqc_PDB/")

#### Plot data loss from CD-Hit ####

read_split_bar_plots <- function(path, source) {
  table_path <- paste0(path, "multiqc_report_data/")
  files <- list.files(table_path, pattern = "multiqc_split_bar_plot_*", full.names = TRUE)
  ids <- sub("^multiqc_split_bar_plot_(.*)\\.txt$", "\\1", basename(files))
  tables <- lapply(files, fread, sep = "\t")
  names(tables) <- ids
  combined <- rbindlist(tables, idcol = "ID")
  combined$Source <- source
  combined
}

cdhit <- rbindlist(list(
  read_split_bar_plots(path_to_results, "non_struc"),
  read_split_bar_plots(path_to_struc_results, "structural"),
  read_split_bar_plots(path_to_pdb_results, "PDB")
))

cdhit$Dataset <- sapply(str_split(cdhit$ID, "-", n = 2), `[`, 1)
cdhit$Category <- sapply(str_split(cdhit$ID, "-", n = 2), `[`, 2)
cdhit[, Category := factor(Category, levels = c("regular", "ILP-split", "ILP-split-and-sampling", "Leakage-split", "ILP-hcn"))]
cdhit <- cdhit[!(Dataset == "PDB_dimers" & Source == "structural")]
cdhit <- cdhit[, -c("ID", "Source")]
cdhit <- cdhit[!Category %in% c("ILP-split-and-sampling", "Leakage-split")]

cdhit_plt <- melt(cdhit, id.vars = c("Sample", "Dataset", "Category"))
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


rename_map_sample <- c(
  "1_train" = "Train",
  "2_val" = "Val", 
  "3_test" = "Test",
  "6_discarded" = "Discarded"
)
cdhit_plt[, Sample := rename_map_sample[Sample]]
cdhit_plt[, Sample := factor(Sample, levels = c("Train", "Val", "Test", "Discarded"))]
cdhit_plt <- cdhit_plt[!Category %in% c("ILP-split-and-sampling", "Leakage-split")]
cdhit_plt[Category == "regular", Category := "Similarity-reduced"]
cdhit_plt[Category == "ILP-split", Category := "ILP split"]

ggplot(cdhit_plt, aes(x = value, y = Sample, fill = variable))+
  geom_col()+
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
ggsave("figures/cdhit_loss.pdf", width=12, height=15)
