# Load necessary libraries
library(VennDiagram)
library(readr)
library(grid)
library(openxlsx)
library(dplyr)

# Directory containing the DESeq2 output files for Anatomic Features
#directory <- "/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/top_proviruses/SHS_AS/FL"
directory = 'C:/Users/worst/iCloudDrive/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/top_proviruses/SHS_AS/FL'


# List all CSV files in the directory
files <- list.files(directory, pattern = "*.csv", full.names = TRUE)

# Function to extract anatomical features from filenames
extract_features <- function(filename) {
  # Extract the part of the filename after the last underscore and before "_FL.csv"
  base <- basename(filename)
  features <- sub(".*_top_proviruses_", "", base)
  features <- sub("_FL.csv", "", features)
  
  # Split features by "_or_" and handle "CTTRUE" case
  feature_list <- unlist(strsplit(features, "_or_"))
  
  # Append "_base" to "CTTRUE"
  feature_list <- sub("CTTRUE", "CT_base", feature_list)
  
  # Remove "TRUE" from feature names
  feature_list <- sub("TRUE", "", feature_list)
  return(feature_list)
}

# Read each file into a list of data frames, naming them by their features
data_list <- lapply(files, function(file) {
  features <- extract_features(file)
  data <- read_csv(file)
  data$features <- list(features)
  return(data)
})

# Combine all data frames into one, with a column for features
all_data <- bind_rows(data_list)

# Function to get HERVs with a specific sign for a given individual feature
get_hervs_individual <- function(data, feature, sign) {
  data %>%
    dplyr::filter(sapply(features, function(f) identical(f, feature))) %>%
    dplyr::filter(sign(log2FoldChange) == sign) %>%
    pull(1) %>%
    unique()
}

# Function to get HERVs with a specific sign for a given pair or three-way feature
get_hervs_combined <- function(data, feature_combination, sign) {
  data %>%
    dplyr::filter(sapply(features, function(f) setequal(f, feature_combination))) %>%
    dplyr::filter(sign(log2FoldChange) == sign) %>%
    pull(1) %>%
    unique()
}

# Define the anatomical features
features <- c("CT_base", "CTmvp", "CTpan", "ITcx", "LEcx")

# Function to assign HERVs to feature intersections
assign_hervs_to_features <- function(data, sign) {
  herv_assignments <- list()
  
  # Track presence in individual, pair, and three-way feature files
  individual_presence <- list()
  pair_presence <- list()
  three_presence <- list()
  
  # Check individual feature files
  for (feature in features) {
    hervs <- get_hervs_individual(data, feature, sign)
    for (herv in hervs) {
      individual_presence[[herv]] <- c(individual_presence[[herv]], feature)
    }
  }
  
  # Check pair feature files
  pair_features <- combn(features, 2, simplify = FALSE)
  for (pair in pair_features) {
    hervs <- get_hervs_combined(data, pair, sign)
    for (herv in hervs) {
      pair_presence[[herv]] <- c(pair_presence[[herv]], pair)
    }
  }
  
  # Check three-way feature files
  three_features <- combn(features, 3, simplify = FALSE)
  for (three in three_features) {
    hervs <- get_hervs_combined(data, three, sign)
    for (herv in hervs) {
      three_presence[[herv]] <- c(three_presence[[herv]], three)
    }
  }
  
  # Assign HERVs based on presence
  for (herv in unique(c(names(individual_presence), names(pair_presence), names(three_presence)))) {
    if (!is.null(individual_presence[[herv]]) && length(individual_presence[[herv]]) == 1) {
      # Assign to the single individual feature
      herv_assignments[[herv]] <- individual_presence[[herv]]
    } else if (!is.null(individual_presence[[herv]])) {
      # Assign to the intersection of individual features
      herv_assignments[[herv]] <- individual_presence[[herv]]
    } else if (!is.null(pair_presence[[herv]]) && length(pair_presence[[herv]]) == 1) {
      # Assign to the single pair feature
      herv_assignments[[herv]] <- unique(unlist(pair_presence[[herv]]))
    } else if (!is.null(pair_presence[[herv]])) {
      # Assign to the intersection of pair features
      pairs <- unique(unlist(pair_presence[[herv]]))
      herv_assignments[[herv]] <- pairs
    } else if (!is.null(three_presence[[herv]]) && length(three_presence[[herv]]) == 1) {
      # Assign to the single three-way feature
      herv_assignments[[herv]] <- unique(unlist(three_presence[[herv]]))
    } else if (!is.null(three_presence[[herv]])) {
      # Assign to the intersection of three-way features
      threes <- unique(unlist(three_presence[[herv]]))
      herv_assignments[[herv]] <- threes
    }
  }
  
  return(herv_assignments)
}

# Assign HERVs to features for positive and negative log2FC
positive_assignments <- assign_hervs_to_features(all_data, 1)
negative_assignments <- assign_hervs_to_features(all_data, -1)

# Create sets for overexpressed and underexpressed HERVs
overexpressed_sets <- lapply(features, function(feature) {
  names(positive_assignments)[sapply(positive_assignments, function(x) feature %in% x)]
})
names(overexpressed_sets) <- features

underexpressed_sets <- lapply(features, function(feature) {
  names(negative_assignments)[sapply(negative_assignments, function(x) feature %in% x)]
})
names(underexpressed_sets) <- features

# Update feature names for plot labels
feature_labels <- c("CT", "MVP", "PAN", "IT", "LE")

# Generate Venn diagram for overexpressed HERVs
venn.plot.over <- venn.diagram(
  x = overexpressed_sets,
  category.names = feature_labels,
  filename = NULL,
  output = TRUE,
  cex = 1.5,  # Increase text size for the numbers
  cat.cex = 1.5,  # Increase text size for the category names
  fill = c("red", "green", "blue", "yellow", "purple"),
  alpha = 0.5,
  label.col = "black",
  print.mode = "raw",  # Show raw counts only
  cat.pos = c(0, 0, 0, 0, 360),  # Adjust positions of category labels
  cat.dist = c(0.05, 0.05, 0.05, 0.05, 0.05)  # Adjust distance of category labels
)

# Draw the Venn diagram and set the font to Helvetica
grid.draw(venn.plot.over)

grob_list <- grid.ls(print = FALSE)

for (grob_name in grob_list$name) {
  if (grepl("text", grob_name)) {
    grid.gedit(grob_name, gp = gpar(fontfamily = "Helvetica"))
  }
}

# Generate Venn diagram for underexpressed HERVs
venn.plot.under <- venn.diagram(
  x = underexpressed_sets,
  category.names = feature_labels,
  filename = NULL,
  output = TRUE,
  cex = 1.5,  # Increase text size for the numbers
  cat.cex = 1.5,  # Increase text size for the category names
  fill = c("red", "green", "blue", "yellow", "purple"),
  alpha = 0.5,
  label.col = "black",
  print.mode = "raw",  # Show raw counts only
  cat.pos = c(0, 0, 0, 0, 360),  # Adjust positions of category labels
  cat.dist = c(0.05, 0.05, 0.05, 0.05, 0.05)  # Adjust distance of category labels
)

# Draw the Venn diagram and set the font to Helvetica
grid.draw(venn.plot.under)

# List all grobs to find the correct paths
grob_list <- grid.ls(print = FALSE)

# Loop through all grobs and change font to Helvetica for text elements
for (grob_name in grob_list$name) {
  if (grepl("text", grob_name)) {
    grid.gedit(grob_name, gp = gpar(fontfamily = "Helvetica"))
  }
}

# Plot the Venn diagrams sequentially in the RStudio plotting pane
grid.newpage()
grid.draw(venn.plot.over)

grid.newpage()
grid.draw(venn.plot.under)


##Create workbooks##
# Function to create a workbook with sheets for each intersection
create_workbook <- function(assignments, filename, all_data, sign) {
  wb <- createWorkbook()
  
  # Filter all_data for the relevant sign
  filtered_data <- all_data %>% dplyr::filter(sign(log2FoldChange) == sign)
  
  # Add a sheet for each feature
  for (feature in features) {
    addWorksheet(wb, feature)
    feature_data <- names(assignments)[sapply(assignments, function(x) identical(x, feature))]
    
    # Extract log2FoldChange and padj for each HERV in this feature
    feature_df <- filtered_data %>%
      dplyr::filter(.[[1]] %in% feature_data) %>%
      select(HERV = 1, log2FoldChange, padj, features)
    
    writeData(wb, feature, feature_df)
  }
  
  # Add sheets for pairwise intersections
  pair_features <- combn(features, 2, simplify = FALSE)
  for (pair in pair_features) {
    intersection_name <- paste(pair, collapse = "_and_")
    addWorksheet(wb, intersection_name)
    intersection_data <- names(assignments)[sapply(assignments, function(x) setequal(x, pair))]
    
    # Extract log2FoldChange and padj for each HERV in this intersection
    intersection_df <- filtered_data %>%
      dplyr::filter(.[[1]] %in% intersection_data) %>%
      select(HERV = 1, log2FoldChange, padj, features)
    
    writeData(wb, intersection_name, intersection_df)
  }
  
  # Add sheets for three-way intersections
  three_features <- combn(features, 3, simplify = FALSE)
  for (three in three_features) {
    intersection_name <- paste(three, collapse = "_and_")
    addWorksheet(wb, intersection_name)
    intersection_data <- names(assignments)[sapply(assignments, function(x) setequal(x, three))]
    
    # Extract log2FoldChange and padj for each HERV in this intersection
    intersection_df <- filtered_data %>%
      dplyr::filter(.[[1]] %in% intersection_data) %>%
      select(HERV = 1, log2FoldChange, padj, features)
    
    writeData(wb, intersection_name, intersection_df)
  }
  
  # Save the workbook
  saveWorkbook(wb, filename, overwrite = TRUE)
}

# Create workbooks for overexpressed and underexpressed HERVs
#create_workbook(positive_assignments, "/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/top_proviruses/SHS_AS/overexpressed_HERVs.xlsx", all_data, 1)
#create_workbook(negative_assignments, "/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/top_proviruses/SHS_AS/underexpressed_HERVs.xlsx", all_data, -1)