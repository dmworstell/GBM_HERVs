# Load the necessary package
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Define the file path
base_dir    = '/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis'

#file_path_DEprvGeneOLs <- file.path(base_dir,"top_proviruses","output_overlapping_pvs.csv")
file_path_DEprvGeneOLs <- "/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs/overlaps_DEprvs_genes.csv"
file_path_AllDEprvs <- "/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/top_proviruses/090623_results_GBM_vs_Controls_top_proviruses_withSense.csv"
#file_path_AllDEprvs <- "/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/top_proviruses/2024-10-22_top_proviruses_AllHervs.csv"
# Read the CSV files into dataframes
DEPRVGeneOLs <- read.csv(file_path_DEprvGeneOLs,header=TRUE)
#names(DEPRVGeneOLs) <- "Provirus"
DEPrvs <- read.csv(file_path_AllDEprvs)
colnames(DEPrvs)[1] <- "Provirus"

# Use the anti_join function to get the rows in DEPrvs that are not in DEPRVGeneOLs
DEPrvsNoGeneOL <- anti_join(DEPrvs, DEPRVGeneOLs, by = "Provirus")

# Create a dataframe with rows where "exp_type" is "underexpressed"
DEPRVGeneOLs_underexpressed <- subset(DEPRVGeneOLs, exp_type == "underexpressed")

# Create a dataframe with rows where "exp_type" is "overexpressed"
DEPRVGeneOLs_overexpressed <- subset(DEPRVGeneOLs, exp_type == "overexpressed")

# Define a function to process a subset
process_subset <- function(df) {
  df$Overlap_Type <- paste(unique(df$transcript_type), collapse = "_AND_")
  df <- df[1, ]
  return(df)
}

# For the underexpressed dataframe
df_list <- split(DEPRVGeneOLs_underexpressed, DEPRVGeneOLs_underexpressed$Provirus)
DEPRVGeneOLs_underexpressed_grouped <- do.call(rbind, lapply(df_list, process_subset))

# For the overexpressed dataframe
df_list <- split(DEPRVGeneOLs_overexpressed, DEPRVGeneOLs_overexpressed$Provirus)
DEPRVGeneOLs_overexpressed_grouped <- do.call(rbind, lapply(df_list, process_subset))

#Adjust labels#

# Define a function to replace specific strings
replace_strings <- function(x) {
  x <- sub("lncRNA_AND_protein_coding", "protein_coding", x)
#  x <- sub("lncRNA_AND_protein_coding", "protein_coding_AND_lncRNA", x)
  x <- sub("protein_coding_AND_lncRNA", "protein_coding", x)
  return(x)
}

# For the DEPRVGeneOLs_underexpressed_grouped dataframe
DEPRVGeneOLs_underexpressed_grouped$Overlap_Type <- replace_strings(DEPRVGeneOLs_underexpressed_grouped$Overlap_Type)

# For the DEPRVGeneOLs_overexpressed_grouped dataframe
DEPRVGeneOLs_overexpressed_grouped$Overlap_Type <- replace_strings(DEPRVGeneOLs_overexpressed_grouped$Overlap_Type)

DEPrvNoGeneOL_OE <- subset(DEPrvsNoGeneOL,log2FoldChange > 0)
DEPrvNoGeneOL_OE$Overlap_Type <- "Intergenic"
DEPrvNoGeneOL_UE <- subset(DEPrvsNoGeneOL,log2FoldChange < 0)
DEPrvNoGeneOL_UE$Overlap_Type <- "Intergenic"

## Add in the proviruses without Gene OL
# Select only the columns you care about from each dataframe
DEPrvNoGeneOL_OE_selected <- subset(DEPrvNoGeneOL_OE, select = c(Provirus, Overlap_Type))
DEPRVGeneOLs_overexpressed_grouped_selected <- subset(DEPRVGeneOLs_overexpressed_grouped, select = c(Provirus, Overlap_Type))

# Use rbind to combine the dataframes
combined_df_OE <- rbind(DEPrvNoGeneOL_OE_selected, DEPRVGeneOLs_overexpressed_grouped_selected)


## Add in the proviruses without Gene OL
# Select only the columns you care about from each dataframe
DEPrvNoGeneOL_UE_selected <- subset(DEPrvNoGeneOL_OE, select = c(Provirus, Overlap_Type))
DEPRVGeneOLs_underexpressed_grouped_selected <- subset(DEPRVGeneOLs_underexpressed_grouped, select = c(Provirus, Overlap_Type))

# Use rbind to combine the dataframes
combined_df_UE <- rbind(DEPrvNoGeneOL_UE_selected, DEPRVGeneOLs_underexpressed_grouped_selected)

##Generate Pie Charts##
library(RColorBrewer)

# Define a colorblind-friendly palette
color_palette <- brewer.pal(n = length(unique(combined_df_UE$Overlap_Type)), name = "Set2")

# Function to replace underscores with spaces and adjust specific labels
replace_labels <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("protein coding", "mRNA", x)
  return(x)
}

###########DEPrvs_UE#########
# Calculate the percentages
overlap_counts <- table(combined_df_UE$Overlap_Type)
overlap_percentages <- overlap_counts / sum(overlap_counts)

# Create a data frame for the percentages
percentage_df <- data.frame(Overlap_Type = names(overlap_percentages), 
                            percentage = as.vector(overlap_percentages))

# Create the pie chart
ggplot(percentage_df, aes(x = "", y = percentage, fill = Overlap_Type)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  theme_void() +
  labs(title = "DEPrvs_UE Overlap Types", fill = "Overlap Type") +
  scale_fill_manual(values = color_palette, labels = replace_labels(names(overlap_percentages))) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

###########DEPrvs_OE#########
# Calculate the percentages
overlap_counts <- table(combined_df_OE$Overlap_Type)
overlap_percentages <- overlap_counts / sum(overlap_counts)

# Create a data frame for the percentages
percentage_df <- data.frame(Overlap_Type = names(overlap_percentages), 
                            percentage = as.vector(overlap_percentages))

# Create the pie chart
ggplot(percentage_df, aes(x = "", y = percentage, fill = Overlap_Type)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  theme_void() +
  labs(title = "DEPrvs_OE Overlap Types", fill = "Overlap Type") +
  scale_fill_manual(values = color_palette, labels = replace_labels(names(overlap_percentages))) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))