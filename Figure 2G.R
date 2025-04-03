###############
#Takes in combined counts data from telescope
#Takes in a file with metadata including variables to compare
#Outputs Differential Expression matrices

#Daniel Murimi-Worstell
#March 3, 2022

#Edited:
#5/15/24
###############

#Load packages
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(apeglm))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(stringr))


#Set up directory paths
Genome="hg38"
base_dir    = '/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis'
results_dir = file.path(base_dir, 'results')

##### Read in counts #####
raw_count_file_dir  = file.path(results_dir,'telescope', 'counts_matrices', 'raw',Genome,'All_Proviruses')
control_counts = read.csv(paste(raw_count_file_dir,'Control_Telescope_output_hg38_p14_RepeatMasker_HS_LS_Reports_LS.csv',sep='/'))
gbm_counts = read.csv(paste(raw_count_file_dir,'Glioblastoma_Telescope_output_hg38_p14_RepeatMasker_HS_LS_Reports_LS.csv',sep='/'))

#Combine counts files
counts <- cbind(gbm_counts,control_counts[,-1])

#Convert counts to matrix
countmat = round(data.matrix(counts[,-1]))
rownames(countmat) = counts$X
rm(counts)

#Read in column data
control_coldata = read.csv(paste(results_dir,'runTable_Control.csv',sep='/'))
gbm_coldata = read.csv(paste(results_dir,'runTable_GBM.csv',sep='/'))

#Remove empty rows
gbm_coldata <- gbm_coldata[!apply(gbm_coldata == "", 1, all), ]
rownames(gbm_coldata) <- NULL

#Combine and reset indexes
coldata <- rbind(gbm_coldata, control_coldata)
rownames(coldata) <- NULL

#Remove rows that weren't run through telescope
coldata <- coldata[coldata$run %in% colnames(countmat), ] 
rownames(coldata) <- NULL

#Removing ALL SHS runs
just_using_runs <- coldata$run[!(coldata$Dataset == 'SHS')]
coldata<- subset(coldata,Dataset!='SHS')

countmat <- countmat[, just_using_runs]

## Read in the names of proviruses that overlap genes so they can be filtered out ##
prv_names <- read.csv(file.path(base_dir,"top_proviruses","output_overlapping_pvs.csv"),header=FALSE)
prv_names <- unlist(prv_names)
## subset the counts matrix accordingly ##
countmat_sub <- countmat[rownames(countmat)[!(rownames(countmat) %in% prv_names)],]

####Decide here whether to use the filtered matrix or just regular matrix
countmat <- countmat_sub
###

#####DESEQ2#####

dds <- DESeqDataSetFromMatrix(countData = countmat,
                              colData   = coldata,
                              design    = ~ disease_state )
#Filter genes
keep = which(rowMeans(counts(dds)) > 2)
dds <- dds[keep,]

#Go to the correct directory
setwd('/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/figures')

#run DEseq2
dds <- DESeq(dds)

#Look at results
resultsNames(dds)
results = lfcShrink(dds, coef = "disease_state_GBM_vs_Control", type = "apeglm")

results = results[order(results$padj),]

#Get dimensions of results
dim(results)

#Check head of results
head(results)

#Filter these genes
top_genes = subset(results, (abs(results$log2FoldChange) > 1.5) & 
                     (-log10(results$padj) > -log10(0.05)))

OE_genes = as.data.frame(subset(results, (results$log2FoldChange > 1.5) & 
                     (-log10(results$padj) > -log10(0.05))))

UE_genes = as.data.frame(subset(results, (results$log2FoldChange < -1.5) & 
                    (-log10(results$padj) > -log10(0.05))))
  
#Write them out
top_genes = as.data.frame(top_genes)

###################################
### NEW METHOD OF VISUALIZATION ###
###################################

###Select exonized intergenic proviruses
# Define the path to your file
file_path <- "/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/top_proviruses/hg38/New_Transcripts/TCGA/intergenic/output_proviruses.tsv"

# Read the file into a character vector
lines <- readLines(file_path)

OE_top_genes <- top_genes[top_genes$log2FoldChange > 1.5,]
top_prvs <- OE_top_genes[rownames(OE_top_genes) %in% lines,]

# Calculate the number of samples in each condition
n_control <- sum(coldata$disease_state == "Control")
n_gbm <- sum(coldata$disease_state == "GBM")

# Calculate the harmonic mean of the sample sizes
harmonic_mean <- 2 * n_control * n_gbm / (n_control + n_gbm)

# Calculate the standard deviations using the harmonic mean
top_prvs$lfcSD <- top_prvs$lfcSE * sqrt(harmonic_mean)

plot_data <- data.frame(
  provirus = rownames(top_prvs),
  log2_FC = as.numeric(top_prvs$log2FoldChange),
  log2_FC_SE = as.numeric(top_prvs$lfcSE),
  log2_FC_SD = as.numeric(top_prvs$lfcSD)
)

# Display the results
print(top_genes)

library(ggplot2)

# Assuming 'plot_data' is your data frame and it has columns 'provirus', 'log2_FC', and 'log2_FC_SD'

# Create a new column to specify the color group
plot_data$color_group <- ifelse(plot_data$provirus == "2q12.1_chr2:103477983-103481010_MER57-int", "MER57-int 2q12.1", "other")

# Reorder the provirus factor based on log2_FC from low to high
plot_data$provirus <- reorder(plot_data$provirus, plot_data$log2_FC)

# Plot
ggplot(plot_data, aes(x = provirus, y = log2_FC)) +
  geom_point(aes(color = color_group), size = 3) +
  geom_errorbar(aes(ymin = log2_FC - log2_FC_SD, ymax = log2_FC + log2_FC_SD), width = 0.2) +
  scale_color_manual(values = c("MER57-int 2q12.1" = "purple", "other" = "red"), 
                     labels = c("MER57-int 2q12.1" = "MER57-int 2q12.1", "other" = "Other")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Exonized HERVs in Novel Transcripts",
    x = NULL,
    y = expression(Log[2] * "(Fold Change)")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Hide x-axis text
    axis.text.y = element_text(size = 12),  # Increase y-axis text size
    axis.ticks.y = element_line(color = "black"),  # Add y-axis tick marks
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Optional: Remove panel background
    plot.margin = margin(1, 1, 1, 1, "cm"),  # Add margins around the plot
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center and bold title
    panel.border = element_rect(color = "black", fill = NA)  # Add axis lines
  ) +
  guides(color = guide_legend(title = NULL))  # Add legend with no title

#Write out sheet for supporting data
library(openxlsx)

# Define the file path
file_path <- '/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Manuscripts/2024/Supplemental/Supporting_Data.xlsx'

# Define the sheet name
sheet_name <- "Figure 2B"

# Check if the file exists
if (file.exists(file_path)) {
  # Load the existing workbook
  wb <- loadWorkbook(file_path)
} else {
  # Create a new workbook
  wb <- createWorkbook()
}

# Add a new sheet with the specified name
addWorksheet(wb, sheet_name)

# Write the data to the new sheet
writeData(wb, sheet_name, plot_data)

# Save the workbook
saveWorkbook(wb, file_path, overwrite = TRUE)