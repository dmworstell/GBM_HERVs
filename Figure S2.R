###############

#Daniel Murimi-Worstell
#December 6 2023

###############

#Load packages
suppressPackageStartupMessages(library(apeglm))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(PCAtools))


#Set up directory paths
Genome="hg38"
#base_dir    = '/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis'
base_dir = 'C:/Users/worst/iCloudDrive/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis'

results_dir = file.path(base_dir, 'results')
raw_count_file_dir  = file.path(results_dir,'telescope', 'counts_matrices', 'raw',Genome,'All_Proviruses')

### LOADING PREVIOUSLY MADE COUNTS MATRIX ###
countmat <- read.csv(file.path(results_dir,file='telescope','counts_matrices','raw','hg38','All_Proviruses','Glioblastoma_Telescope_output_hg38_p14_RepeatMasker_HS_LS_Reports_LS.csv'),row.names=1)
counts <- as.data.frame(countmat)
#######

#Read in column data
coldata = read.csv(paste(results_dir,'runTable_GBM.csv',sep='/'))

#Remove empty rows
coldata <- coldata[!apply(coldata == "", 1, all), ]
rownames(coldata) <- NULL

#Remove rows that weren't run through telescope
coldata <- coldata[coldata$run %in% names(counts), ] 
rownames(coldata) <- NULL

#Only TCGA Runs
just_using_runs <- coldata$run[coldata$Dataset == 'TCGA']
coldata<- subset(coldata,Dataset=='TCGA')
countmat <- countmat[, just_using_runs]

## Read in the names of proviruses that overlap genes so they can be filtered out ##
prv_names <- read.csv(file.path(base_dir,"top_proviruses","output_overlapping_pvs.csv"),header=FALSE)
prv_names <- unlist(prv_names)
## subset the counts matrix accordingly ##
countmat_sub <- countmat[rownames(countmat)[!(rownames(countmat) %in% prv_names)],]

####Decide here whether to use the filtered matrix or just regular matrix
countmat <- countmat_sub
###


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
                              design    = ~ subtype )

#Filter genes
keep = which(rowMeans(counts(dds)) > 2)
dds <- dds[keep,]

#Go to the correct directory
setwd('C:/Users/worst/iCloudDrive/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/figures')

#run DEseq2
dds <- DESeq(dds)

#Variance Stabilize
vsdds = varianceStabilizingTransformation(dds)
vs_counts = assays(vsdds)[[1]]
md <- coldata
rownames(md) <- colnames(vs_counts)
p <- pca(vs_counts, metadata = md)

####
##EIGENCORPLOT##

elbow <- findElbowPoint(p$variance)

to_factor <- c("recurrence","subtype","gender","IDH1_status","MGMT_PCR","PTEN",
               'EGFR','EGFR_vIII','Radiation','Chemotherapy')

make_numeric <- c('OS','PFS','Percent.aneuploidy','TERT.expression.log2',
                  'ABSOLUTE.purity','ABSOLUTE.ploidy','Age_Diagnosis')
p_factored <- p
p_factored[["metadata"]][to_factor] <- lapply(p[["metadata"]][to_factor], factor)
p_factored[["metadata"]][make_numeric] <- lapply(p[["metadata"]][make_numeric], as.numeric)

# Rename the columns in the metadata
colnames(p_factored[["metadata"]]) <- gsub("disease_state", "Disease State", colnames(p_factored[["metadata"]]))
colnames(p_factored[["metadata"]]) <- gsub("IDH1_status", "IDH1 Mutation", colnames(p_factored[["metadata"]]))
colnames(p_factored[["metadata"]]) <- gsub("MGMT_PCR", "MGMT Promoter Methylation", colnames(p_factored[["metadata"]]))
colnames(p_factored[["metadata"]]) <- gsub("recurrence", "Recurrence", colnames(p_factored[["metadata"]]))
colnames(p_factored[["metadata"]]) <- gsub("gender", "Gender", colnames(p_factored[["metadata"]]))
colnames(p_factored[["metadata"]]) <- gsub("subtype", "Subtype", colnames(p_factored[["metadata"]]))

colnames(p_factored[["metadata"]]) <- gsub("ABSOLUTE.ploidy", "Ploidy", colnames(p_factored[["metadata"]]))
colnames(p_factored[["metadata"]]) <- gsub("ABSOLUTE.purity", "Purity", colnames(p_factored[["metadata"]]))
colnames(p_factored[["metadata"]]) <- gsub("TERT.expression.log2", "TERT Expression", colnames(p_factored[["metadata"]]))
colnames(p_factored[["metadata"]]) <- gsub("Age_Diagnosis", "Age at Diagnosis", colnames(p_factored[["metadata"]]))
colnames(p_factored[["metadata"]]) <- gsub("EGFR_vIII", "EGFR vIII", colnames(p_factored[["metadata"]]))

# Define the metavars with the new names
mvs <- c('Subtype', 'IDH1 Mutation', 'Recurrence', 'Gender', 
         'MGMT Promoter Methylation', 'PTEN', 'EGFR', 
         'Radiation', 'Chemotherapy', 'OS', 'PFS',
         'Ploidy','Purity','TERT Expression','Age at Diagnosis','EGFR vIII')

eigencorplot(p_factored,
             components = getComponents(p, 1:20),
             metavars = mvs)