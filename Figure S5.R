suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

base_dir <- "/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/results/bootstrapping_results"

# Function to calculate combined median ratios
calculate_combined_median_ratios <- function(de_data, all_data) {
  # Combine categories
  de_combined <- data.frame(
    Exon = rowSums(de_data[, c("+ | First", "- | First", "+ | Middle", "- | Middle", "+ | Last", "- | Last", "+ | Single", "- | Single")], na.rm = TRUE),
    Intron = rowSums(de_data[, c("+ | Intron", "- | Intron")], na.rm = TRUE)
  )
  
  all_combined <- data.frame(
    Exon = rowSums(all_data[, c("+ | First", "- | First", "+ | Middle", "- | Middle", "+ | Last", "- | Last", "+ | Single", "- | Single")], na.rm = TRUE),
    Intron = rowSums(all_data[, c("+ | Intron", "- | Intron")], na.rm = TRUE)
  )
  
  # Calculate medians for combined categories
  de_medians <- apply(de_combined, 2, median)
  all_medians <- apply(all_combined, 2, median)
  
  # Calculate ratios of medians
  ratios <- de_medians / all_medians
  return(ratios)
}

# Calculate ratios for each threshold
thresholds <- c("100", "80", "60", "40")
lncRNA_ratios <- list()
mRNA_ratios <- list()

for (threshold in thresholds) {
  # Load data for the current threshold
  lncRNA_DE_prvs <- read.csv(file.path(base_dir, paste0("lncRNA_proviruses_1000_bootstraps_7146_MT_", threshold, "percentExon_DEHervs.csv")), row.names = 1, check.names=FALSE)
  lncRNA_All_Proviruses <- read.csv(file.path(base_dir, paste0("lncRNA_proviruses_1000_bootstraps_7146_MT_", threshold, "percentExon_AllHervs.csv")), row.names = 1, check.names=FALSE)
  mRNA_DE_prvs <- read.csv(file.path(base_dir, paste0("mRNA_proviruses_1000_bootstraps_7146_MT_", threshold, "percentExon_DEHervs.csv")), row.names = 1, check.names=FALSE)
  mRNA_All_Proviruses <- read.csv(file.path(base_dir, paste0("mRNA_proviruses_1000_bootstraps_7146_MT_", threshold, "percentExon_AllHervs.csv")), row.names = 1, check.names=FALSE)
  
  # Calculate combined median ratios
  lncRNA_ratios[[threshold]] <- calculate_combined_median_ratios(lncRNA_DE_prvs, lncRNA_All_Proviruses)
  mRNA_ratios[[threshold]] <- calculate_combined_median_ratios(mRNA_DE_prvs, mRNA_All_Proviruses)
}

# Combine ratios into a data frame for plotting
plot_data <- data.frame()
for (threshold in thresholds) {
  lncRNA_df <- data.frame(
    Threshold = threshold,
    Category = names(lncRNA_ratios[[threshold]]),
    Ratio = lncRNA_ratios[[threshold]],
    Type = "lncRNA"
  )
  mRNA_df <- data.frame(
    Threshold = threshold,
    Category = names(mRNA_ratios[[threshold]]),
    Ratio = mRNA_ratios[[threshold]],
    Type = "mRNA"
  )
  plot_data <- rbind(plot_data, lncRNA_df, mRNA_df)
}

# Debugging: Print the plot_data to check its structure
print(head(plot_data))

# Set the order of the Threshold factor
plot_data$Threshold <- factor(plot_data$Threshold, levels = c("40", "60", "80", "100"))

# Plot the data
ggplot(plot_data, aes(x = Category, y = Ratio, fill = Threshold)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ Type) +
  theme(
    axis.text = element_text(size = 12),  # Increase font size for axis text
    axis.title = element_text(size = 14), # Increase font size for axis titles
    plot.title = element_text(size = 16, hjust = 0.5), # Increase font size and center the title
    strip.text = element_text(size = 14)  # Increase font size for facet labels
  ) +
  labs(title = "Relative Enrichment of DE to All HERVs",
       y = "Relative Enrichment (DE : All)",
       x = "Overlap Category",
       fill = "Threshold")