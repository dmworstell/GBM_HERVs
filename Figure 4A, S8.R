###############
#Takes in combined counts data from telescope
#Takes in a file with metadata including variables to compare
#Outputs Differential Expression matrices

#Daniel Murimi-Worstell
#September 26 2023

#Edited:
#10/22/24
###############

#Load packages
suppressPackageStartupMessages(library(pheatmap))
#suppressPackageStartupMessages(library(genefilter))
suppressPackageStartupMessages(library(apeglm))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(PCAtools))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyverse))


#Set up directory paths
Genome="hg38"
#base_dir    = '/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis'
base_dir = 'C:/Users/worst/iCloudDrive/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis'
results_dir = file.path(base_dir, 'results')
raw_count_file_dir  = file.path(results_dir,'telescope', 'counts_matrices', 'raw',Genome,'All_Proviruses')

### LOADING PREVIOUSLY MADE COUNTS MATRIX ###
countmat <- read.csv(file.path(results_dir,file='telescope','counts_matrices','raw','hg38','All_Proviruses','Glioblastoma_Telescope_output_hg38_p14_RepeatMasker_HS_LS_Reports_LS.csv'),row.names=1)
counts <- as.data.frame(countmat)
#######

#Read in column data
coldata = read.csv(paste(results_dir,'runTable_SHS_Only.csv',sep='/'))

#Remove empty rows
coldata <- coldata[!apply(coldata == "", 1, all), ]
rownames(coldata) <- NULL

#Remove rows that weren't run through telescope
coldata <- coldata[coldata$run %in% names(counts), ] 
rownames(coldata) <- NULL

#Removing TCGA runs
#just_using_runs <- coldata$run[!(coldata$study_name == 'TCGA')]
just_using_runs <- coldata$run
coldata<- subset(coldata,study_name!='TCGA')
countmat <- countmat[, just_using_runs]

#Removing Stem Cell runs
just_using_runs <- coldata$run[!(coldata$study_name == 'CSC')]
coldata<- subset(coldata,study_name!='CSC')
countmat <- countmat[, just_using_runs]

#Removing CTmvp samples
#just_using_runs <- coldata$run[!(coldata$structure_parent_acronym == 'CTmvp')]
#coldata<- subset(coldata,structure_parent_acronym!='CTmvp')
#countmat <- countmat[, just_using_runs]

#Removing LEcx samples
#just_using_runs <- coldata$run[!(coldata$structure_parent_acronym == 'LEcx')]
#coldata<- subset(coldata,structure_parent_acronym!='LEcx')
#countmat <- countmat[, just_using_runs]

#Removing CTpan samples
#just_using_runs <- coldata$run[!(coldata$structure_parent_acronym == 'CTpan')]
#coldata<- subset(coldata,structure_parent_acronym!='CTpan')
#countmat <- countmat[, just_using_runs]

# Create a new coldata with TRUE/FALSE values for each unique value in the structure_parent_acronym column
coldata_binary <- coldata %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = structure_parent_acronym, values_from = value, values_fill = FALSE)
coldata_binary$structure_parent_acronym <- coldata$structure_parent_acronym

# Create a list of all unique pairs of structure_parent_acronym values
spa_pairs <- combn(unique(coldata$structure_parent_acronym), 2, simplify = FALSE)

# Create new columns for each pair
for(i in seq_along(spa_pairs)) {
  pair <- spa_pairs[[i]]
  colname <- paste(pair, collapse = "_or_")
  coldata_binary[[colname]] <- coldata_binary[[pair[1]]] | coldata_binary[[pair[2]]]
}

# Create a list of all unique trios of structure_parent_acronym values
spa_trios <- combn(unique(coldata$structure_parent_acronym), 3, simplify = FALSE)

# Create new columns for each trio
for(i in seq_along(spa_trios)) {
  trio <- spa_trios[[i]]
  colname <- paste(trio, collapse = "_or_")
  coldata_binary[[colname]] <- coldata_binary[[trio[1]]] | coldata_binary[[trio[2]]] | coldata_binary[[trio[3]]]
}

## Read in the names of proviruses that overlap genes so they can be filtered out ##
prv_names <- read.csv(file.path(base_dir,"top_proviruses","output_overlapping_pvs.csv"),header=FALSE)
prv_names <- unlist(prv_names)
## subset the counts matrix accordingly ##
countmat_sub <- countmat[rownames(countmat)[!(rownames(countmat) %in% prv_names)],]

####Decide here whether to use the filtered matrix or just regular matrix
countmat <- countmat_sub

###
# Make life easier by filtering out proviruses with median counts under 2 here ##

#Filter proviruses
keep = which(rowMeans(countmat) > 2)
countmat <- countmat[keep,]


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


#dds <- DESeqDataSetFromMatrix(countData = countmat,
#                             colData   = coldata_binary,
#                             design   = ~ structure_parent_acronym + subtype)

#Go to the correct directory
setwd('C:/Users/worst/iCloudDrive/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/figures')

toLoop <-colnames(coldata_binary[21:45])
toLoop <- toLoop[!grepl("CSC_Other",toLoop)]

for (var in toLoop) {
  design_formula <- as.formula(paste("~ ", var, "+ subtype"))
  dds <- DESeqDataSetFromMatrix(countData = countmat,
                                colData   = coldata_binary,
                                design = design_formula)
#  dds <- DESeq(dds)
  
#dds <- DESeqDataSetFromMatrix(countData = countmat,
#                             colData   = coldata_binary,
#                             design = ~ CT + subtype )
    

#dds <- DESeqDataSetFromMatrix(countData = countmat,
#                             colData   = coldata_binary,
#                             design   = ~ CT )


#run DEseq2
dds <- DESeq(dds)

#Look at results
resultsNames(dds)

coefToUse = resultsNames(dds)[2]
results = lfcShrink(dds, coef = coefToUse , type = "apeglm")

results = results[order(results$padj),]


top_proviruses = subset(results, (abs(results$log2FoldChange) > 1.5) & 
                     (-log10(results$padj) > -log10(0.05)))

OE_proviruses = as.data.frame(subset(results, (results$log2FoldChange > 1.5) & 
                     (-log10(results$padj) > -log10(0.05))))

UE_proviruses = as.data.frame(subset(results, (results$log2FoldChange < -1.5) & 
                    (-log10(results$padj) > -log10(0.05))))
  
#Write them out
top_proviruses = as.data.frame(top_proviruses)

#Write results
write.csv(top_proviruses, file = file.path(base_dir,'top_proviruses/SHS_AS/',paste0(Sys.Date(),'_SHS_AS_top_proviruses_',coefToUse, '.csv')),  
          row.names = TRUE, quote = FALSE)

write.csv(results, file = file.path(base_dir,'results/SHS_AS/',paste0(Sys.Date(),'_SHS_AS_results_',coefToUse,'.csv')),  
                    row.names = TRUE, quote = FALSE)
}

#Plot out specific provirus
#vsdds = vst(dds)
vsdds = varianceStabilizingTransformation(dds)
vs_counts = assays(vsdds)[[1]]
# 
### Make Heatmap ###
##### hg38 proviruses ######
meta = colData(dds)[,c("structure_parent_acronym","subtype")]

meta = as.data.frame(meta)

rownames(meta) <- coldata$run

#SORT meta
meta_sorted <- meta[order(meta[,1],meta[,2], decreasing=FALSE),]

provirusmat <- vs_counts

#pvsofinterest <- pvsofinterest[pvsofinterest %in% rownames(countmat)]

#OVEREXPRESSED
#CT OR LEcx, LTR45B
#CT OR ITcx, LTR88a
#                   CTmvp                           CTpan                       
#pvsofinterest <- c("LTR59","LTR18C","MER52D","LTR3B_v","MER54B","HERVIP10F-int","LTR45B",
#                   "LTR50","MER61-int","LTR48B","MER31B","LTR3","MER57C2","MER51-int","LTR88a")
#                   #LEcx AND ITcx                   #LEcx NOT ITcx               From previous figure, ITcx, CT
 
#FROM PRIOR FIGURE    
#pvsofinterest <- c("LTR3B_v","HERVIP10F-int","LTR102_Mam","MER57C2","MER51-int",
#                  "LTR88a",  "LTR48B","MER61-int","LTR50","LTR3","LTR10B1","MER41D",
#                  "LTR33","LTR1A1","HERVL-int","MER21-int","MER4B-int")


#HIGH EXPRESSION
pvsofinterest <- c( "LTR3B_v","HERVIP10F-int","LTR102_Mam", #PAN
                   "LTR52","LTR4","LTR45B", #CT_base
                   "MER66D","LTR9C","LTR27C", #ITcx
                   "LTR48B","MER61-int","LTR50", #LEcx
                   "LTR33","LTR1A1","MER21-int") #MVP
#LOW EXPRESSION

# pvsofinterest <- c( "MamGypLTR2b","MER39B","HERVL-int" #PAN
#                     "LTR3","MamGypsy2-I","MER52D", #CT_base
#                     "MER4B-int","MamGypLTR2c","HERVK3-int",#ITcx
#                     "MLT1E3","MER57-int","MER73", #LEcx
#                     "LTR48B","MER61-int","LTR75_1") #MVP

selectedmat <- provirusmat[pvsofinterest,]

#SORT selectedmat
selectedmat_sorted <- selectedmat[,rownames(meta_sorted)]

pheatmap(selectedmat_sorted, cluster_rows = FALSE, 
         show_rownames = TRUE, show_colnames=FALSE, cluster_cols = FALSE, annotation_col=meta_sorted)




# Create a named vector for renaming the categories
rename_map <- c(
  "CT" = "CT",
  "CTmvp" = "MVP",
  "LEcx" = "LE",
  "ITcx" = "IT",
  "CTpan" = "PAN"
)

# Ensure meta_sorted is a data frame with row names matching the heatmap columns
annotation_col <- data.frame(structure_parent_acronym = meta_sorted$structure_parent_acronym)
rownames(annotation_col) <- rownames(meta_sorted)

# Apply the renaming to the structure_parent_acronym column
annotation_col$structure_parent_acronym <- rename_map[annotation_col$structure_parent_acronym]

# Rename the column in the annotation data frame to "Anatomic Feature"
colnames(annotation_col) <- "Anatomic Feature"

# Define the desired order of the Anatomic Features
desired_order <- c("PAN", "CT", "IT", "LE", "MVP")

# Rename the structure_parent_acronym column in meta_sorted

meta_sorted$structure_parent_acronym <- rename_map[meta_sorted$structure_parent_acronym]

# Reorder meta_sorted based on the desired order
meta_sorted <- meta_sorted[order(match(meta_sorted$structure_parent_acronym, desired_order)), ]

# Reorder the selectedmat_sorted matrix columns according to the new order
selectedmat_sorted <- selectedmat_sorted[, rownames(meta_sorted)]

# Create the annotation_col data frame
annotation_col <- data.frame(Anatomic_Feature = meta_sorted$structure_parent_acronym)
rownames(annotation_col) <- rownames(meta_sorted)

annotation_col$Anatomic_Feature <- factor(annotation_col$Anatomic_Feature, levels = desired_order)

# Create the heatmap with the specified customizations
pheatmap(selectedmat_sorted, 
         cluster_rows = FALSE, 
         show_rownames = TRUE, 
         show_colnames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = annotation_col,
         annotation_names_col = FALSE,  # Remove the label above the heatmap
         fontsize = 12)  # Increase the size of all text, including the legend

md <- coldata
rownames(md) <- coldata$run

md <- md %>%
  mutate(structure_parent_acronym = recode(structure_parent_acronym,
                                           'CT' = 'CT',
                                           'CTmvp' = 'MVP',
                                           'CTpan' = 'PAN',
                                           'ITcx' = 'IT',
                                           'LEcx' = 'LE'))


p <- pca(vs_counts, metadata = md)

# Define a named vector of colors for your categories with hex code placeholders
category_colors <- c(
  "CT" = "#ff9695",             
  "IT" = "#fffd96",        
  "LE" = "#d298f7",                
  "MVP" = "#9ffb95", 
  "PAN" = "#9397ff"             
)

# Use the colkey argument to specify the colors
bp <- biplot(p,
       showLoadings = T,
       lengthLoadingsArrowsFactor = 2,
       showLoadingsNames = F,
       #sizeLoadingsNames = 0,
       x = 'PC1', y = 'PC2',
       xlim = c(-12, 16),
       ylim = c(-11, 7),
       axisLabSize = 20,  # increase axis label text size
       lab = p$run,
       colby = 'structure_parent_acronym',
       colkey = category_colors,  # Use the named vector here
       hline = 0, vline = 0,
       legendPosition = 'none',
       title = "Clades by Anatomic Feature")  # Add the title here


png("Figure4A.png", width = 5, height = 5, units = "in", res = 300, type = "cairo")
print(bp)
dev.off()