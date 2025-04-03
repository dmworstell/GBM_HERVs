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

control_counts = read.csv(paste(raw_count_file_dir, 'Control_Telescope_output_hg38_p14_RepeatMasker_HS_LS_Reports_LS.csv',sep='/'))
gbm_counts = read.csv(paste(raw_count_file_dir,'Glioblastoma_Telescope_output_hg38_p14_RepeatMasker_HS_LS_Reports_LS.csv',sep='/'))

#Combine counts files
counts <- cbind(gbm_counts,control_counts[,-1])

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
coldata <- coldata[coldata$run %in% names(counts), ] 
rownames(coldata) <- NULL

#Convert counts to matrix
countmat = round(data.matrix(counts[,-1]))
rownames(countmat) = counts$X
rm(counts)

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
setwd('/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/figures')

#run DEseq2
dds <- DESeq(dds)

#Look at results
resultsNames(dds)
results = lfcShrink(dds, coef = "disease_state_GBM_vs_Control", type = "apeglm")

results = results[order(results$padj),]
results

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

#Write results
#write.csv(top_genes, file = file.path(base_dir,'top_proviruses',paste0(Sys.Date(),'_top_proviruses_AllHervs.csv')),  
#          row.names = TRUE, quote = FALSE)

#write.csv(results, file = file.path(base_dir,'results',paste0(Sys.Date(),'_results_GBM_vs_Controls_AllHervs.csv')),  
#          row.names = TRUE, quote = FALSE)

# Generate volcano plot
de_genes = subset(results, abs(log2FoldChange) > 1.5 & padj < 0.05)

# Assign colors based on log2FoldChange
# High fold change (positive) will be red, low fold change (negative) will be blue
de_color = ifelse(results$log2FoldChange > 1.5, "red", 
                  ifelse(results$log2FoldChange < -1.5, "blue", "black"))

# Plot with custom axis labels
plot(results$log2FoldChange, -log10(results$padj), pch = 20, col = de_color,
     xlab = "Log2(Fold Change)", ylab = "-Log10(Adjusted p-values)")


vsdds = varianceStabilizingTransformation(dds)
vs_counts = assays(vsdds)[[1]]


md <- coldata
rownames(md) <- colnames(vs_counts)

# Create a new column in the metadata that combines IDH1_status and disease_state
md$combined_status <- paste(md$IDH1_status, md$disease_state, sep = "_")

# Ensure the new column is a factor
md$combined_status <- factor(md$combined_status,
                             levels = c("Mutant_GBM", "Wildtype_Control", "Wildtype_GBM"),
                             labels = c("IDH1mt GBM", "Controls", "IDH1wt GBM"))

p <- pca(vs_counts, metadata = md)

library(RColorBrewer)

# Define a color palette with transparency
# Set transparency using the alpha channel (values between 0 and 1, where 0 is fully transparent and 1 is fully opaque)
palette_colors <- c(
  adjustcolor(brewer.pal(n = 3, name = "Set2")[3], alpha.f = 1),   # IDH1mt GBM (fully opaque)
  adjustcolor(brewer.pal(n = 3, name = "Set2")[1], alpha.f = 0.3), # Controls (semi-transparent)
  adjustcolor(brewer.pal(n = 3, name = "Set2")[2], alpha.f = 0.3)  # IDH1wt GBM (semi-transparent)
)

# Assuming 'p' is your data frame and 'combined_status' is a factor
#p$combined_status_f <- factor(p$combined_status,
#                            levels = c("Mutant_GBM", "Wildtype_Control", "Wildtype_GBM"),
#                            labels = c("IDH1mt", "Controls", "IDH1wt GBM"))


##BOTH##
# biplot(p,
#         showLoadings = FALSE,
#        #lengthLoadingsArrowsFactor = 1.5,
#        #sizeLoadingsNames = 4,
#        x='PC1',y='PC2',
#        lab = p$run,
#        colby = c('disease_state'),
#        colkey = palette_colors,
#        hline = 0, vline = 0,
#        legendPosition = 'right')


##To show IDHmut status##

# Generate the biplot using the combined status for coloring
biplot(p,
       showLoadings = FALSE,
       x = 'PC1', y = 'PC2',
       lab = p$run,
       colby = 'combined_status',  # Use the updated combined column
       colkey = palette_colors,
       hline = 0, vline = 0,
       legendPosition = 'right',
       title = "Clades by IDH1 status")

##EIGENCORPLOT##
elbow <- findElbowPoint(p$variance)

to_factor <- c("disease_state","recurrence","subtype","gender","IDH1_status","MGMT_PCR","PTEN",
               'EGFR','Radiation','Chemotherapy')

make_numeric <- c('OS','PFS')

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

# Define the metavars with the new names
mvs <- c('Disease State', 'Subtype', 'IDH1 Mutation', 'Recurrence', 'Gender', 
         'MGMT Promoter Methylation', 'PTEN', 'EGFR', 
         'Radiation', 'Chemotherapy', 'OS', 'PFS')

# Create the eigencorplot with updated metadata
eigencorplot(p_factored,
             components = getComponents(p, 1:20),
             metavars = mvs)

