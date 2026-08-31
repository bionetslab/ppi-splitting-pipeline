library(data.table)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(fmsb)

path_to_results <- "~/PythonProjects/ppi-splitting-pipeline/results_server/multiqc/"
path_to_struc_results <- "~/PythonProjects/ppi-splitting-pipeline/results_server/multiqc_structural/"


test_results_balanced <- fread(paste0(path_to_results, "multiqc_report_data/multiqc_classifier_metrics_test_balanced_table.txt"), sep="\t")
test_results_balanced_struc <- fread(paste0(path_to_struc_results, "multiqc_report_data/multiqc_classifier_metrics_test_balanced_table.txt"), sep="\t")
test_results_balanced <- rbind(test_results_balanced, test_results_balanced_struc)
test_results_balanced$Dataset <- sapply(stringr::str_split(test_results_balanced$Sample, "-", n = 2),`[`, 1)
test_results_balanced$Category <- sapply(stringr::str_split(test_results_balanced$Sample, "-", n = 2),`[`, 2)
test_results_balanced[, Category := factor(Category, levels = c("Leakage-split", "regular", "ILP-split", "ILP-split-and-sampling"))]
test_results_balanced[, Dataset := factor(Dataset, levels = c("PDB_dimers", "Pinder", "Intact", "Hippie", "String", "String900", "StringDB", "StringExp", "StringText"))]

test_results_real <- fread(paste0(path_to_results, "multiqc_report_data/multiqc_classifier_metrics_test_realistic_table.txt"), sep="\t")
test_results_real_struc <- fread(paste0(path_to_struc_results, "multiqc_report_data/multiqc_classifier_metrics_test_realistic_table.txt"), sep="\t")
test_results_real <- rbind(test_results_real, test_results_real_struc)
test_results_real$Dataset <- sapply(stringr::str_split(test_results_real$Sample, "-", n = 2),`[`, 1)
test_results_real$Category <- sapply(stringr::str_split(test_results_real$Sample, "-", n = 2),`[`, 2)
test_results_real[, Category := factor(Category, levels = c("Leakage-split", "regular", "ILP-split", "ILP-split-and-sampling"))]
test_results_real[, Dataset := factor(Dataset, levels = c("PDB_dimers", "Pinder", "Intact", "Hippie", "String", "String900", "StringDB", "StringExp", "StringText"))]

test_results_balanced$TestSet <- "balanced"
test_results_real$TestSet <- "realistic"

only_leakage <- rbind(test_results_balanced, test_results_real)
only_leakage <- only_leakage[Category == "Leakage-split"]
only_leakage <- only_leakage[, -c("Sample", "ID", "Category")]
only_leakage <- only_leakage[order(AUROC)]
only_leakage <- only_leakage[order(Dataset, TestSet)]

leakage_mat <- as.matrix(only_leakage[, c("AUROC", "AUPRC", "F1", "Precision", "Recall", "Accuracy", "MCC")])
rownames(leakage_mat) <- paste(only_leakage$Dataset, only_leakage$TestSet)

compute_number_colors <- function(mat) {   
  breaks <- pheatmap:::generate_breaks(mat, 100)    
  color <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  # extract the hex code that would be assigned to each of your numbers and convert it into rgb 
  rgb_colors <- col2rgb(pheatmap:::scale_colours(as.matrix(mat), col = color, breaks = breaks, na_col = "#DDDDDD"))
  # calculate the value of the function
  luminance <- rgb_colors * c(0.299, 0.587, 0.114)
  luminance <- luminance['red', ]+luminance['green', ] + luminance['blue', ]
  # apply the threshold - I went for a softer black and white to reduce eye strain
  number_color <- ifelse(luminance < 125, "grey90", "grey30")
  return(number_color) 
}

number_color = compute_number_colors(leakage_mat)

pheatmap(leakage_mat,
         cluster_rows = F, 
         cluster_cols = F,
         display_numbers = T, 
         gaps_row = c(2, 4, 6, 8, 10, 12, 14, 16),
         number_color = number_color, 
         breaks = pheatmap:::generate_breaks(leakage_mat, 100),
         legend = F)

path_to_tables <- "multiqc/multiqc_report_data/"

all_tables <- lapply(list.files(paste0(
  path_to_results, path_to_tables), 
  pattern="multiqc_split_bar_plot_*", full.names = T), 
  function(x){fread(x, sep="\t")})

col_names <- list.files(paste0(path_to_results, path_to_tables), pattern="multiqc_split_bar_plot_*")
col_names <- tstrsplit(unname(sapply(col_names, function(x){tstrsplit(x, '_', keep=5)[[1]]})), '\\.', keep=1)[[1]]
names(all_tables) <- col_names
sizes <- rbindlist(all_tables, idcol = "ID")
sizes$Dataset <- sapply(stringr::str_split(sizes$ID, "-", n = 2),`[`, 1)
sizes$Category <- sapply(stringr::str_split(sizes$ID, "-", n = 2),`[`, 2)
sizes[, Sample := tstrsplit(Sample, "_", keep=2)]

sizes <- melt(sizes, id.vars = c("ID", "Sample", "Dataset", "Category"))
sizes[, Sample := factor(Sample, levels=c('discarded', 'train', 'val', 'test'))]
sizes[, Category := factor(Category, levels = c('regular', 'ILP-split', 'ILP-split-and-sampling', 'Leakage-split'))]

ggplot(sizes, aes(x = Category, y = value, fill = variable))+
  facet_grid(Dataset~Sample, scales = 'free')+
  geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))
