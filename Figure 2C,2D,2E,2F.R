# Load the necessary package
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)

# Define the file path
file_path <- "/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs/overlaps_DEprvs_genes.csv"

# Read the CSV file into a dataframe
DEPrvGeneOLs <- read_csv(file_path)

# Create a dataframe with rows where "exp_type" is "underexpressed"
DEPrvGeneOLs_underexpressed <- subset(DEPrvGeneOLs, exp_type == "underexpressed")

#Get list of UE HERVs in final exon
FinalExonDEPrvs_UE <- unique(select(DEPrvGeneOLs_underexpressed[DEPrvGeneOLs_underexpressed$is_last_exon == TRUE,],c('Provirus','sense','Strand','gene_name','transcript_type','relative_orientation')))

# Create a dataframe with rows where "exp_type" is "overexpressed"
DEPrvGeneOLs_overexpressed <- subset(DEPrvGeneOLs, exp_type == "overexpressed")

#Get list of OE HERVs in final exon
FinalExonDEPrvs_OE <- unique(select(DEPrvGeneOLs_overexpressed[DEPrvGeneOLs_overexpressed$is_last_exon == TRUE,],c('Provirus','sense','Strand','gene_name','transcript_type','relative_orientation')))

# Function to add new columns
add_columns <- function(df) {
  df %>%
    mutate(
      Chromosome = str_extract(Provirus, "(?<=_)[^:]+"),
      Start = str_extract(Provirus, "(?<=:)[^\\-]+"),
      End = str_extract(Provirus, "(?<=-)[^_]+")
    )
}

# Apply the function to both dataframes
FinalExonDEPrvs_OE <- add_columns(FinalExonDEPrvs_OE)
FinalExonDEPrvs_UE <- add_columns(FinalExonDEPrvs_UE)

#Fix alt assembly provirus chromosome names
FinalExonDEPrvs_OE <- FinalExonDEPrvs_OE %>% mutate(
  Chromosome = sub("v",".",str_extract(Chromosome, "[^_]+")))

FinalExonDEPrvs_UE <- FinalExonDEPrvs_UE %>% mutate(
  Chromosome = sub("v",".",str_extract(Chromosome, "[^_]+")))

# Write out the resulting dataframes of final-exon-overlapping HERVs
#write.csv(FinalExonDEPrvs_OE,"/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs/FinalExonDEPrvsOverExpressed.csv",row.names=FALSE)
#write.csv(FinalExonDEPrvs_UE,"/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs/FinalExonDEPrvsUnderExpressed.csv",row.names=FALSE)

# Define a function to process a subset
process_subset <- function(df) {
  df$Overlap_Type <- paste(unique(df$exon_pos), collapse = "_AND_")
  df <- df[1, ]
  return(df)
}

# For the underexpressed dataframe
df_list <- split(DEPrvGeneOLs_underexpressed, DEPrvGeneOLs_underexpressed$Provirus)
DEPrvGeneOLs_underexpressed_grouped <- do.call(rbind, lapply(df_list, process_subset))

# For the overexpressed dataframe
df_list <- split(DEPrvGeneOLs_overexpressed, DEPrvGeneOLs_overexpressed$Provirus)
DEPrvGeneOLs_overexpressed_grouped <- do.call(rbind, lapply(df_list, process_subset))

##Split into PC and LNCRNA
# For the underexpressed dataframe
PrvGnOLs_UE_PC <- subset(DEPrvGeneOLs_underexpressed_grouped, transcript_type == "protein_coding")
PrvGnOLs_UE_lncRNA <- subset(DEPrvGeneOLs_underexpressed_grouped, transcript_type == "lncRNA")

# For the overexpressed dataframe
PrvGnOLs_OE_PC <- subset(DEPrvGeneOLs_overexpressed_grouped, transcript_type == "protein_coding")
PrvGnOLs_OE_lncRNA <- subset(DEPrvGeneOLs_overexpressed_grouped, transcript_type == "lncRNA")


#Adjust labels#

# Define a function to replace specific strings
replace_strings <- function(x) {
  x <- sub("intron_AND_", "", x)
  x <- sub("middle_exon_AND_first_exon", "first_exon", x)
  x <- sub("first_exon_AND_middle_exon", "first_exon", x)
  x <- sub("first_exon_AND_single_exon", "first_exon", x)
  x <- sub("middle_exon_AND_last_exon", "last_exon", x)
  x <- sub("last_exon_AND_middle_exon", "last_exon", x)
  x <- sub("last_exon_AND_single_exon", "last_exon", x)
  x <- sub("first_exon_AND_last_exon_AND_middle_exon", "first_exon_AND_last_exon", x)
  x <- sub("first_exon_AND_last_exon_AND_single_exon", "first_exon_AND_last_exon", x)
  x <- sub("single_exon_AND_last_exon", "last_exon", x)
  x <- sub("single_exon_AND_first_exon", "first_exon", x)
  x <- sub("single_exon_AND_middle_exon", "middle_exon", x)
  x <- sub("single_exon_AND_intron", "single_exon", x)
  x <- sub("last_exon_AND_first_exon", "first_exon_AND_last_exon", x)
  
  return(x)
}

# For the PrvGnOLs_UE_PC dataframe
PrvGnOLs_UE_PC$Overlap_Type <- replace_strings(PrvGnOLs_UE_PC$Overlap_Type)

# For the PrvGnOLs_UE_lncRNA dataframe
PrvGnOLs_UE_lncRNA$Overlap_Type <- replace_strings(PrvGnOLs_UE_lncRNA$Overlap_Type)

# For the PrvGnOLs_OE_PC dataframe
PrvGnOLs_OE_PC$Overlap_Type <- replace_strings(PrvGnOLs_OE_PC$Overlap_Type)

# For the PrvGnOLs_OE_lncRNA dataframe
PrvGnOLs_OE_lncRNA$Overlap_Type <- replace_strings(PrvGnOLs_OE_lncRNA$Overlap_Type)

##Generate Pie Charts##

###########PrvGnOLs_UE_PC#########
# Calculate the percentages
overlap_counts <- table(PrvGnOLs_UE_PC$Overlap_Type)
overlap_percentages <- overlap_counts / sum(overlap_counts)

# Create a data frame for the percentages
percentage_df <- data.frame(Overlap_Type = names(overlap_percentages), 
                            percentage = as.vector(overlap_percentages))

# Create the pie chart
ggplot(percentage_df, aes(x = "", y = percentage, fill = Overlap_Type)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  #geom_text_repel(aes(label = scales::percent(percentage, accuracy = 0.01)), position = position_stack(vjust = 0.5)) +
  theme_void() +
  labs(title = "PrvGnOLs_UE_PC Overlap Types")

###########PrvGnOLs_UE_lncRNA#########
# Calculate the percentages
overlap_counts <- table(PrvGnOLs_UE_lncRNA$Overlap_Type)
overlap_percentages <- overlap_counts / sum(overlap_counts)

# Create a data frame for the percentages
percentage_df <- data.frame(Overlap_Type = names(overlap_percentages), 
                            percentage = as.vector(overlap_percentages))

# Create the pie chart
ggplot(percentage_df, aes(x = "", y = percentage, fill = Overlap_Type)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  #geom_text_repel(aes(label = scales::percent(percentage, accuracy = 0.01)), position = position_stack(vjust = 0.5)) +
  theme_void() +
  labs(title = "PrvGnOLs_UE_lncRNA Overlap Types")

###########PrvGnOLs_OE_PC#########
# Calculate the percentages
overlap_counts <- table(PrvGnOLs_OE_PC$Overlap_Type)
overlap_percentages <- overlap_counts / sum(overlap_counts)

# Create a data frame for the percentages
percentage_df <- data.frame(Overlap_Type = names(overlap_percentages), 
                            percentage = as.vector(overlap_percentages))

# Create the pie chart
ggplot(percentage_df, aes(x = "", y = percentage, fill = Overlap_Type)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  #geom_text_repel(aes(label = scales::percent(round(percentage, 2))), position = position_stack(vjust = 0.5)) +
  theme_void() +
  labs(title = "PrvGnOLs_OE_PC Overlap Types")


###########PrvGnOLs_OE_lncRNA#########
# Calculate the percentages
overlap_counts <- table(PrvGnOLs_OE_lncRNA$Overlap_Type)
overlap_percentages <- overlap_counts / sum(overlap_counts)

# Create a data frame for the percentages
percentage_df <- data.frame(Overlap_Type = names(overlap_percentages), 
                            percentage = as.vector(overlap_percentages))

# Create the pie chart
ggplot(percentage_df, aes(x = "", y = percentage, fill = Overlap_Type)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  #geom_text_repel(aes(label = scales::percent(percentage, accuracy = 0.01)), position = position_stack(vjust = 0.5)) +
  theme_void() +
  labs(title = "PrvGnOLs_OE_lncRNA Overlap Types")