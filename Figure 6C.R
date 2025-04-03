library(tidyverse)
#library(correlationAnalyzeR)
library(DGCA)
suppressPackageStartupMessages(library(apeglm))
suppressPackageStartupMessages(library(DESeq2))
library(ggplot2, quietly = TRUE)
library(GOstats, quietly = TRUE)
library(HGNChelper, quietly = TRUE)
library(org.Hs.eg.db, quietly = TRUE)
library(plotrix, quietly = TRUE)

#Load in data
#base_dir    = '/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis'
base_dir = 'C:/Users/worst/iCloudDrive/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis'

results_dir = file.path(base_dir, 'results')
Genome="hg38"
raw_count_file_dir  = file.path(results_dir,'telescope', 'counts_matrices', 'raw',Genome,'HML2_Genes')
counts = read.csv(paste(raw_count_file_dir,'Glioblastoma_Telescope_output_HML2_GENES.csv',sep='/'))

###Get the HERV sheet to combine###
raw_count_file_dir_HERV  = file.path(results_dir,'telescope', 'counts_matrices', 'raw',Genome,'All_Proviruses')
GBM_counts_HERV = read.csv(paste(raw_count_file_dir_HERV,'Glioblastoma_Telescope_output_hg38_p14_RepeatMasker_HS_LS_Reports_LS.csv',sep='/'))

geneOfInterest = "2q12.1_chr2:103477983-103481010_MER57-int"

GBM_counts_HERV_chk <- subset(GBM_counts_HERV,X==geneOfInterest)
counts <- rbind(counts,GBM_counts_HERV_chk)
#counts <- rbind(counts, GBM_counts_HERV)

rm(GBM_counts_HERV)


#Convert counts to matrix
countmat = round(data.matrix(counts[,-1]))
rownames(countmat) = counts$X
rm(counts)

#Read in column data
coldata = read.csv(paste(results_dir,'runTable_GBM.csv',sep='/'))

#Remove empty rows
coldata <- coldata[!apply(coldata == "", 1, all), ]
rownames(coldata) <- NULL

#Remove rows that weren't run through telescope
coldata <- coldata[coldata$run %in% colnames(countmat), ] 
rownames(coldata) <- NULL

#Just SHS_AS runs
just_using_runs <- coldata$run[(coldata$study_name == 'AnatomicStructures')]
coldata<- subset(coldata,study_name=='AnatomicStructures')
countmat <- countmat[, just_using_runs]

#Just SHS runs
#just_using_runs <- coldata$run[(coldata$Dataset == 'SHS')]
#coldata<- subset(coldata,Dataset=='SHS')
#countmat <- countmat[, just_using_runs]


#Filter
keep = which(rowMeans(countmat) > 2)
countmat <- countmat[keep,]

### Analyze correlations

#dds <- DESeqDataSetFromMatrix(countData = countmat,
#                              colData   = coldata,
#                              design    = ~ subtype + structure_parent_acronym)

dds <- DESeqDataSetFromMatrix(countData = countmat,
                              colData   = coldata,
                              design    = ~ subtype + structure_parent_acronym)


#dds <- DESeqDataSetFromMatrix(countData = countmat,
#                              colData   = coldata,
#                              design    = ~ subtype)


dds <- DESeq(dds)

vsdds = vst(dds)
vs_counts = assays(vsdds)[[1]]

vs_counts_cv_filtered = filterGenes(vs_counts, 
                                   filterTypes = c("central", "dispersion"), filterDispersionType = "cv", 
                                   filterDispersionPercentile = 0.1, sequential = TRUE)


vs_counts <- vs_counts_cv_filtered
rm(vs_counts_cv_filtered)

###MAKE SELECTIONS HERE
Category <- "subtype"
Grp1 <- "Classical"
Grp2 <- "Mesenchymal"
#Category <- "structure_parent_acronym"
#Grp1 <- "CTpan"
#Grp2 <- "LEcx"
######

n_Grp1_samples = length(coldata[[Category]][coldata[[Category]] == Grp1])
n_Grp2_samples = length(coldata[[Category]][coldata[[Category]] == Grp2]) 
sample_type = c(rep(Grp1, n_Grp1_samples),rep(Grp2, n_Grp2_samples))
design_mat = makeDesign(sample_type)


### REORDER AND FILTER###
# Filter coldata to keep only rows where the specified category equals Grp1 or Grp2
coldata_filtered <- coldata[coldata[[Category]] %in% c(Grp1, Grp2), ]

# Reorder the rows of coldata by the specified category
coldata_filtered <- coldata_filtered[order(coldata_filtered[[Category]]), ]

# Reorder the columns of vs_counts to match the new order of the run column in coldata
vs_counts_filtered <- vs_counts[, coldata_filtered$run]

###

ddcor_res_GOI = ddcorAll(inputMat = vs_counts_filtered, corrType = "spearman", design = design_mat,
                         compare = c(Grp1, Grp2), sigThresh=0.05,
                         adjust = "BH", heatmapPlot = FALSE, nPerm = 0, splitSet = geneOfInterest)

ddcorGO_res = ddcorGO(ddcor_res_GOI, universe = rownames(vs_counts_filtered), 
                      gene_ontology = "all", HGNC_clean = TRUE, HGNC_switch = TRUE, annotation = "org.Hs.eg.db", calculateVariance = TRUE, adjusted = TRUE)

ddcorGO_BP_gain_of_cor_df = ddcorGO_res[["enrichment_significant_gain_of_correlation_genes"]][["BP"]]
ddcorGO_BP_loss_of_cor_df = ddcorGO_res[["enrichment_significant_loss_of_correlation_genes"]][["BP"]]

### PLOTTING ###

plotCors(inputMat = vs_counts_filtered, design = design_mat, compare = c(Grp1, Grp2), geneA = geneOfInterest, geneB = "KIF14")

sig_gain_df_list <- ddcorGO_res[["enrichment_significant_gain_of_correlation_genes"]]
sig_loss_df_list <- ddcorGO_res[["enrichment_significant_loss_of_correlation_genes"]]

# Extract the "BP" data frames
gain_bp_df <- sig_gain_df_list[["BP"]]
loss_bp_df <- sig_loss_df_list[["BP"]]

# Sort the gain "BP" data frame by OddsRatio in descending order
gain_bp_sorted <- gain_bp_df[order(-gain_bp_df$OddsRatio), ]

# Sort the loss "BP" data frame by OddsRatio in descending order
loss_bp_sorted <- loss_bp_df[order(-loss_bp_df$OddsRatio), ]

# Replace the original "BP" data frames with the sorted versions
sorted_sig_gain_df_list <- sig_gain_df_list
sorted_sig_loss_df_list <- sig_loss_df_list

sorted_sig_gain_df_list[["BP"]] <- gain_bp_sorted
sorted_sig_loss_df_list[["BP"]] <- loss_bp_sorted

# Plot using the sorted lists
labels = c(
  as.expression(substitute(paste("Classical ", rho ," < "," Mesenchymal ", rho), list(Grp1 = Grp1, Grp2 = Grp2))),
  "GO Term",
  as.expression(substitute(paste("Classical ", rho ," > "," Mesenchymal ", rho), list(Grp1 = Grp1, Grp2 = Grp2)))
)

#plotGOTwoGroups(sorted_sig_gain_df_list, sorted_sig_loss_df_list, labels = labels, GOTermTypes = c("BP"))


 custom_order <- c(
   "regulation of chemotaxis",
   "positive regulation of cell adhesion",
   "negative regulation of immune system process",
   "regulation of leukocyte migration",
   "regulation of cell adhesion",
   "chromosome organization",
   "cell cycle process",
   "cell division",
   "mitotic cell cycle",
   "regulation of cell cycle"
 )
 
plotGOTwoGroupsCustom(dfList1 = sorted_sig_gain_df_list, dfList2 = sorted_sig_loss_df_list, labels = labels, GOTermTypes = c("BP"), customOrder = custom_order)