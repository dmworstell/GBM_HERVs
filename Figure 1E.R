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
base_dir    = '/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis'
results_dir = file.path(base_dir, 'results')

##### Read in counts #####
raw_count_file_dir  = file.path(results_dir,'telescope', 'counts_matrices', 'raw',Genome,'HML2_Genes')
control_counts = read.csv(paste(raw_count_file_dir,'Control_Telescope_output_HML2_GENES.csv',sep='/'))
gbm_counts = read.csv(paste(raw_count_file_dir,'Glioblastoma_Telescope_output_HML2_GENES.csv',sep='/'))

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

#Remove all IDHmut samples
just_using_runs <- coldata$run[coldata$IDH1_status == 'Wildtype']
coldata<- subset(coldata,IDH1_status=='Wildtype')
countmat <- countmat[, just_using_runs]

#####DESEQ2#####

dds <- DESeqDataSetFromMatrix(countData = countmat,
                              colData   = coldata,
                              design    = ~ disease_state )
#Filter genes
keep = which(rowMeans(counts(dds)) > 2)
dds <- dds[keep,]

#Go to the correct directory
setwd('/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/figures')

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

pvsofinterest <- c(rownames(top_genes)[grepl("HML-2",rownames(top_genes))])

HERVK_HML2s <-top_genes[rownames(top_genes) %in% pvsofinterest,]

# Calculate the number of samples in each condition
n_control <- sum(coldata$disease_state == "Control")
n_gbm <- sum(coldata$disease_state == "GBM")

# Calculate the harmonic mean of the sample sizes
harmonic_mean <- 2 * n_control * n_gbm / (n_control + n_gbm)

# Calculate the standard deviations using the harmonic mean
HERVK_HML2s$lfcSD <- HERVK_HML2s$lfcSE * sqrt(harmonic_mean)

plot_data <- data.frame(
  provirus = rownames(HERVK_HML2s),
  log2_FC = as.numeric(HERVK_HML2s$log2FoldChange),
  log2_FC_SE = as.numeric(HERVK_HML2s$lfcSE),
  log2_FC_SD = as.numeric(HERVK_HML2s$lfcSD)
)

#colnames(plot_data) <- c("provirus","log2_FC","log2_FC_SE")

# Create the plot
ggplot(plot_data, aes(x = reorder(provirus, -log2_FC), y = log2_FC)) +
  geom_point(aes(color = log2_FC > 1), size = 3) +
  geom_errorbar(aes(ymin = log2_FC - log2_FC_SD, ymax = log2_FC + log2_FC_SD), width = 0.2) +
  scale_color_manual(values = c("blue", "red")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  labs(title = "Log2(Fold Change), Clade Level",
       x = NULL,
       y = "Median log2(FC)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(
    axis.text.x = element_blank(),  # Hide x-axis text
    axis.text.y = element_text(size = 12),  # Increase y-axis text size
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Optional: Remove panel background
    plot.margin = margin(1, 1, 1, 1, "cm")  # Add margins around the plot
  ) +
  guides(color = "none") +
  scale_x_discrete(limits = function(x) c("", x, "")) +  # Add empty categories on both ends
  coord_cartesian(clip = "off")  # Prevent clipping at plot boundaries