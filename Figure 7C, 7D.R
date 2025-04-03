###############
#Takes in combined counts data from telescope
#Takes in a file with metadata including variables to compare
#Outputs Differential Expression matrices

#Daniel Murimi-Worstell
#March 3, 2022

#Edited:
#2/27/23
###############

#Load packages
suppressPackageStartupMessages(library("AnnotationDbi"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(genefilter))
suppressPackageStartupMessages(library(apeglm))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(plyr))
#suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(PCAtools))
#suppressPackageStartupMessages(library(clusterExperiment))
#suppressPackageStartupMessages(library(rgl))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))




Genome="hg38"
base_dir    = '/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis'
results_dir = file.path(base_dir, 'results')
raw_count_file_dir  = file.path(results_dir,'telescope', 'counts_matrices', 'raw',Genome)

###Get GBM data ##
GBM_counts_genes = read.csv(paste(raw_count_file_dir,"HML2_Genes",'Glioblastoma_Telescope_output_HML2_GENES.csv',sep='/'),row.names=1)
GBM_counts_prvs = read.csv(paste(raw_count_file_dir,"All_Proviruses",'Glioblastoma_Telescope_output_hg38_p14_RepeatMasker_HS_LS_Reports_LS.csv',sep='/'),row.names=1)
GBM_counts_genes <- GBM_counts_genes[,colnames(GBM_counts_genes) %in% colnames(GBM_counts_prvs)]
cm_GBM <- rbind(GBM_counts_genes,GBM_counts_prvs)


##Get control data##
control_counts_genes = read.csv(paste(raw_count_file_dir,"HML2_Genes",'Control_Telescope_output_HML2_GENES.csv',sep='/'),row.names=1)
control_counts_prvs = read.csv(paste(raw_count_file_dir,"All_Proviruses",'Control_Telescope_output_hg38_p14_RepeatMasker_HS_LS_Reports_LS.csv',sep='/'),row.names=1)
control_counts_genes <- control_counts_genes[,colnames(control_counts_genes) %in% colnames(control_counts_prvs)]
cm_cntrl <- rbind(control_counts_genes,control_counts_prvs)

countmat <- cbind(cm_GBM, cm_cntrl)

countmat_prvs <- cbind(GBM_counts_prvs,control_counts_prvs)
countmat_genes <- cbind(GBM_counts_genes,control_counts_genes)

rm(GBM_counts_prvs,control_counts_prvs,GBM_counts_genes,control_counts_genes)
#save(countmat,file="/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/results/telescope/counts_matrices/raw/hg38/Combined/Combined.Rdata")

##Load in previously saved combined counts file ##
#load(file.path(results_dir,file='telescope','counts_matrices','raw','hg38','Combined','Combined.Rdata'))

#Read in column data
control_coldata = read.csv(paste(results_dir,'runTable_Control.csv',sep='/'))
gbm_coldata = read.csv(paste(results_dir,'runTable_GBM.csv',sep='/'))

#Combine and reset indexes
coldata <- rbind(gbm_coldata, control_coldata)
rownames(coldata) <- NULL

#Remove rows that weren't run through telescope
coldata <- coldata[coldata$run %in% colnames(countmat), ] 
rownames(coldata) <- NULL

#Removing ALL SHS runs
just_using_runs <- coldata$run[!(coldata$Dataset == 'SHS')]
coldata<- subset(coldata,Dataset!='SHS')

#Removing GTEx Frontal_Cortex runs
#just_using_runs <- just_using_runs[!(coldata$study_name == 'Frontal_Cortex')]
#coldata<- subset(coldata,study_name != 'Frontal_Cortex')


countmat <- countmat[, just_using_runs]
countmat_genes <- countmat_genes[, just_using_runs]
countmat_prvs <- countmat_prvs[, just_using_runs]

rm(control_coldata)
rm(gbm_coldata)



gene_of_interest <- "CARKD"


## Do prvs first
dds <- DESeqDataSetFromMatrix(countData = countmat_prvs,
                              colData   = coldata,
                              design    = ~ disease_state )

#Filter genes
keep = which(rowMeans(counts(dds)) > 2)
dds <- dds[keep,]

#Run DESeq
dds <- DESeq(dds)


vsdds <- vst(dds)
vs_counts_prvs = assays(vsdds)[[1]]

##Then do genes
dds <- DESeqDataSetFromMatrix(countData = countmat_genes,
                              colData   = coldata,
                              design    = ~ disease_state )

#Filter genes
keep = which(rowMeans(counts(dds)) > 2)
dds <- dds[keep,]

#Run DESeq
dds <- DESeq(dds)

vsdds <- vst(dds)
vs_counts_genes = assays(vsdds)[[1]]

##BOX PLOT to show expression in GBM versus controls for a single gene

prv_of_interest <- "13q34_chr13:110613521-110614348_LTR12F"
gene_of_interest <- "CARKD"

prv_exp <- vs_counts_prvs[prv_of_interest,]
gene_exp <- vs_counts_genes[gene_of_interest,]


library(ggplot2)

# Convert named vectors to data frames
prv_exp_df <- data.frame(sample = names(prv_exp), expression = prv_exp)
gene_exp_df <- data.frame(sample = names(gene_exp), expression = gene_exp)

# Merge with metadata
prv_exp_merged <- merge(prv_exp_df, coldata, by.x = "sample", by.y = "run")
gene_exp_merged <- merge(gene_exp_df, coldata, by.x = "sample", by.y = "run")

# Perform Wilcoxon rank-sum test
prv_exp_pvalue <- wilcox.test(expression ~ disease_state, data = prv_exp_merged)$p.value
gene_exp_pvalue <- wilcox.test(expression ~ disease_state, data = gene_exp_merged)$p.value

# Boxplot for prv_exp
ggplot(prv_exp_merged, aes(x = disease_state, y = expression)) +
  geom_boxplot(width = 0.5) +  # Adjust box width
  labs(title = "LTR12F",
       x = NULL,
       y = "Norm. Trans. Expression") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, color="black",face="bold"),
        axis.title = element_text(size = 16,color="black", face="bold"),
        axis.text = element_text(size = 20,color="black"),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines
# +  # Adjust plot margins
#  annotate("text", x = 1.5, y = max(prv_exp_merged$expression), 
#           label = paste("p-value =", signif(prv_exp_pvalue, digits = 3)), 
#           size = 4, vjust = -1)


library(openxlsx)

# Define the file path
file_path <- '/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Manuscripts/2024/Supplemental/Supporting_Data.xlsx'

# Define the sheet name
sheet_name <- "Figure 7B"

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
writeData(wb, sheet_name, prv_exp_merged)

# Save the workbook
saveWorkbook(wb, file_path, overwrite = TRUE)


# Boxplot for gene_exp
ggplot(gene_exp_merged, aes(x = disease_state, y = expression)) +
  geom_boxplot(width = 0.5) +  # Adjust box width
  labs(title = "CARKD",
       x = NULL,
       y = "Norm. Trans. Expression") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, color="black", face="bold"),
        axis.title = element_text(size = 16,color="black", face="bold"),
        axis.text = element_text(size = 20,color="black"),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines
# +  # Adjust plot margins
#  annotate("text", x = 1.5, y = max(gene_exp_merged$expression), 
#           label = paste("p-value =", signif(gene_exp_pvalue, digits = 3)), 
#           size = 4, vjust = -1)


# Define the sheet name
sheet_name <- "Figure 7D"

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
writeData(wb, sheet_name, gene_exp_merged)

# Save the workbook
saveWorkbook(wb, file_path, overwrite = TRUE)