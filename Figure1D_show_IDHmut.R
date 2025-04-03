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
base_dir    = '/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis'
results_dir = file.path(base_dir, 'results')
raw_count_file_dir  = file.path(results_dir,'telescope', 'counts_matrices', 'raw',Genome,'All_Proviruses')

#Read in column data
control_coldata = read.csv(paste(results_dir,'runTable_Control.csv',sep='/'))
gbm_coldata = read.csv(paste(results_dir,'runTable_GBM.csv',sep='/'))

#Remove empty rows
gbm_coldata <- gbm_coldata[!apply(gbm_coldata == "", 1, all), ]
rownames(gbm_coldata) <- NULL

#Combine and reset indexes
coldata <- rbind(gbm_coldata, control_coldata)
rownames(coldata) <- NULL

control_counts = read.csv(paste(raw_count_file_dir, 'Control_Telescope_output_hg38_p14_RepeatMasker_HS_LS_Reports_LS.csv',sep='/'))
gbm_counts = read.csv(paste(raw_count_file_dir,'Glioblastoma_Telescope_output_hg38_p14_RepeatMasker_HS_LS_Reports_LS.csv',sep='/'))


#Combine counts files
counts <- cbind(gbm_counts,control_counts[,-1])

#Remove rows that weren't run through telescope
coldata <- coldata[coldata$run %in% names(counts), ] 
rownames(coldata) <- NULL

#Convert counts to matrix
countmat = round(data.matrix(counts[,-1]))
rownames(countmat) = counts$X
rm(counts)

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
setwd('/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/figures')

#run DEseq2
dds <- DESeq(dds)

#Look at results
#resultsNames(dds)

#Variance Stabilize
vsdds = varianceStabilizingTransformation(dds)
vs_counts = assays(vsdds)[[1]]

md <- coldata
rownames(md) <- colnames(vs_counts)

# Create a new column in the metadata that combines IDH1_status and disease_state
md$combined_status <- paste(md$IDH1_status, md$subtype, sep = "_")

# Ensure the new column is a factor
md$combined_status <- factor(md$combined_status)

# Rename the column in the metadata
colnames(md)[colnames(md) == "IDH1_status"] <- "IDH1 Status"

# Ensure the new column is a factor
md$`IDH1 Status` <- factor(md$`IDH1 Status`)

#PCA
p <- pca(vs_counts, metadata = md)

####
##EIGENCORPLOT##

elbow <- findElbowPoint(p$variance)

to_factor <- c("recurrence","subtype","gender","IDH1_status","MGMT_PCR","PTEN",
               'EGFR','Radiation','Chemotherapy')

make_numeric <- c('OS','PFS','Percent.aneuploidy','TERT.expression.log2',
                  'ABSOLUTE.purity','ABSOLUTE.ploidy')
p_factored <- p
p_factored[["metadata"]][to_factor] <- lapply(p[["metadata"]][to_factor], factor)
p_factored[["metadata"]][make_numeric] <- lapply(p[["metadata"]][make_numeric], as.numeric)

mvs <- c('subtype','IDH1_status','recurrence','gender','MGMT_PCR','ABSOLUTE.purity',
'PTEN','EGFR','Radiation','Chemotherapy','OS','PFS')

eigencorplot(p_factored,
              components = getComponents(p, 1:20),
              metavars = mvs)

library(RColorBrewer)

# Choose a colorblind-friendly palette
palette_colors <- brewer.pal(n = 8, name = "Set2")[c(3,2 )]

palette_colors <- c(
  adjustcolor(brewer.pal(n = 3, name = "Set2")[3], alpha.f = 1),   # IDH1mt GBM (fully opaque)
  adjustcolor(brewer.pal(n = 3, name = "Set2")[2], alpha.f = 0.3)  # IDH1wt GBM (semi-transparent)
)


##PCA ##
# biplot(p,
#        showLoadings = FALSE,
#        #lengthLoadingsArrowsFactor = 1.5,
#        #sizeLoadingsNames = 4,
#        x='PC1',y='PC2',
#        lab = p$run,
#        colby = c('subtype'),
#        #colkey = palette_colors,
#        hline = 0, vline = 0,
#        legendPosition = 'right')

biplot(p,
      showLoadings = FALSE,
      #lengthLoadingsArrowsFactor = 1.5,
      #sizeLoadingsNames = 4,
      x='PC1',y='PC2',
      lab = p$run,
      colby = c('IDH1 Status'),
      colkey = palette_colors,
      hline = 0, vline = 0,
      legendPosition = 'right',
      title="Clades by IDH1 Status")


