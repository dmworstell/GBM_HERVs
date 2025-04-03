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

#Remove all IDHmut samples
just_using_runs <- coldata$run[coldata$IDH1_status == 'Wildtype']
coldata<- subset(coldata,IDH1_status=='Wildtype')
countmat <- countmat[, just_using_runs]

## Read in the names of proviruses that overlap genes so they can be filtered out ##
prv_names <- read.csv(file.path(base_dir,"top_proviruses","output_overlapping_pvs.csv"),header=FALSE)
prv_names <- unlist(prv_names)
## subset the counts matrix accordingly ##
countmat_sub <- countmat[rownames(countmat)[!(rownames(countmat) %in% prv_names)],]

####Decide here whether to use the filtered matrix or just regular matrix
countmat <- countmat_sub
###

### Now, combining into family-level
pvdf <- as.data.frame(countmat)

# Add rownames as a column
pvdf$provirus <- rownames(pvdf)

# Extract family group from provirus column
pvdf$family_group <- gsub(".*-\\d+_", "", pvdf$provirus)

# Reshape from wide to long format
pvdf_long <- reshape2::melt(pvdf, id.vars = c("provirus", "family_group"))

# Sum counts by family group and sample
pvdf_sum <- aggregate(value ~ family_group + variable, pvdf_long, sum)

# Reshape back to wide format
pvdf_wide <- reshape2::dcast(pvdf_sum, family_group ~ variable, value.var = "value")

# Remove the temporary columns from the original dataframe
pvdf$provirus <- NULL
pvdf$family_group <- NULL

rownames(pvdf_wide) <- pvdf_wide$family_group
pvdf_wide$family_group <- NULL
# Print the result
print(pvdf_wide[1:3,1:3])

countmat_backup <- countmat

## Convert count matrix to the family level counts ##
countmat <- pvdf_wide
##

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

pvsofinterest <- c(rownames(top_genes))

# Calculate the number of samples in each condition
n_control <- sum(coldata$disease_state == "Control")
n_gbm <- sum(coldata$disease_state == "GBM")

# Calculate the harmonic mean of the sample sizes
harmonic_mean <- 2 * n_control * n_gbm / (n_control + n_gbm)

# Calculate the standard deviations using the harmonic mean
top_genes$lfcSD <- top_genes$lfcSE * sqrt(harmonic_mean)

plot_data <- data.frame(
  provirus = pvsofinterest,
  log2_FC = as.numeric(top_genes$log2FoldChange),
  log2_FC_SE = as.numeric(top_genes$lfcSE),
  log2_FC_SD = as.numeric(top_genes$lfcSD)
)

# Create the plot
ggplot(plot_data, aes(x = reorder(provirus, -log2_FC), y = log2_FC)) +
  geom_point(aes(color = log2_FC > 1), size = 3) +
  geom_errorbar(aes(ymin = log2_FC - log2_FC_SD, ymax = log2_FC + log2_FC_SD), width = 0.2) +
  scale_color_manual(values = c("blue", "red")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  labs(title = "Log2(Fold Change), Clade Level",
       x = "Provirus Clade",
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