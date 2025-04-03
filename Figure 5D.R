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

Genome="hg38"
results_dir = file.path(base_dir, 'results')
raw_count_file_dir  = file.path(results_dir,'telescope', 'counts_matrices', 'raw',Genome,'HML2_Genes')

control_counts = read.csv(paste(raw_count_file_dir,'Control_Telescope_output_HML2_GENES.csv',sep='/'))
gbm_counts = read.csv(paste(raw_count_file_dir,'Glioblastoma_Telescope_output_HML2_GENES.csv',sep='/'))

###Get the HERV sheet
raw_count_file_dir_HERV  = file.path(results_dir,'telescope', 'counts_matrices', 'raw',Genome,'All_Proviruses')
control_counts_HERV = read.csv(paste(raw_count_file_dir_HERV,'Control_Telescope_output_hg38_p14_RepeatMasker_HS_LS_Reports_LS.csv',sep='/'))
control_counts <- control_counts[,intersect(names(control_counts), names(control_counts_HERV))]

#Combine counts files
counts <- cbind(gbm_counts,control_counts[,-1])

#counts <- gbm_counts
#Convert counts to matrix
countmat = round(data.matrix(counts[,-1]))
rownames(countmat) = counts$X
rm(counts)

#Read in column data
control_coldata = read.csv(paste(results_dir,'runTable_Control.csv',sep='/'))
gbm_coldata = read.csv(paste(results_dir,'runTable_GBM.csv',sep='/'))
#gbm_coldata = subset(gbm_coldata, select = -c(X))

#Remove empty rows
gbm_coldata <- gbm_coldata[!apply(gbm_coldata == "", 1, all), ]
rownames(gbm_coldata) <- NULL

#Pare down to only some columns
#control_coldata <- control_coldata %>% dplyr::select(c('run','disease_state','study_name','subtype','structure_parent_acronym','chr_1p19q_deletion','IDH1_status','recurrence','gender','MGMT_PCR','Dataset','location',
#                                                       'PTEN','EGFR','EGFR_vIII','T2T_PC1_TCGA_Group','Radiation','Chemotherapy','OS','PFS','Race','Age_Diagnosis','Clusters_2012'))

#gbm_coldata <- gbm_coldata %>% dplyr::select(c('run','disease_state','study_name','subtype','structure_parent_acronym','chr_1p19q_deletion','IDH1_status','recurrence','gender','MGMT_PCR','Dataset','location',
#                                               'PTEN','EGFR','EGFR_vIII','T2T_PC1_TCGA_Group','Radiation','Chemotherapy','OS','PFS','Race','Age_Diagnosis','Clusters_2012'))

#coldata <- gbm_coldata

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

#Removing all control runs
#just_using_runs <- coldata$run[!(coldata$disease_state == 'Control')]
#coldata<- subset(coldata,(disease_state!='Control'))
#countmat <- countmat[, just_using_runs]

#Just TCGA runs
#just_using_runs <- coldata$run[(coldata$study_name == 'TCGA')]
#coldata<- subset(coldata,study_name=='TCGA')
#countmat <- countmat[, just_using_runs]

#Filter
keep = which(rowMeans(countmat) > 2)
countmat <- countmat[keep,]

### Analyze correlations

dds <- DESeqDataSetFromMatrix(countData = countmat,
                              colData   = coldata,
                              design    = ~ disease_state )

dds <- DESeq(dds)

#Filter genes
#keep = which(rowMeans(counts(dds)) > 2)
#dds <- dds[keep,]

#resultsNames(dds)
#results = lfcShrink(dds, coef = "disease_state_GBM_vs_Control", type = "apeglm")
#results = results[order(results$padj),]

vsdds = vst(dds)
vs_counts = assays(vsdds)[[1]]

vs_counts_cv_filtered = filterGenes(vs_counts, 
                                   filterTypes = c("central", "dispersion"), filterDispersionType = "cv", 
                                   filterDispersionPercentile = 0.1, sequential = TRUE)
#nrow(vs_counts_cv_filtered)


vs_counts <- vs_counts_cv_filtered
rm(vs_counts_cv_filtered)
#corrMat <- generateCorrelations(countmat,cores=4)

n_ctrl_samples = 252; n_GBM_samples = 169 
#n_ctrl_samples = 252; n_GBM_samples = 444
sample_type = c(rep("GBM", n_GBM_samples),rep("Control", n_ctrl_samples))
design_mat = makeDesign(sample_type)


geneOfInterest = c("HML-2_1q22")

ddcor_res_GOI = ddcorAll(inputMat = vs_counts, corrType = "spearman", design = design_mat,
                         compare = c("GBM", "Control"), sigThresh=0.05,
                         adjust = "BH", heatmapPlot = FALSE, nPerm = 0, splitSet = geneOfInterest)

ddcorGO_res = ddcorGO(ddcor_res_GOI, universe = rownames(vs_counts), 
                      gene_ontology = "all", HGNC_clean = TRUE, HGNC_switch = TRUE, annotation = "org.Hs.eg.db", calculateVariance = TRUE, adjusted = TRUE)

ddcorGO_BP_gain_of_cor_df = ddcorGO_res[["enrichment_significant_gain_of_correlation_genes"]][["BP"]]
ddcorGO_BP_loss_of_cor_df = ddcorGO_res[["enrichment_significant_loss_of_correlation_genes"]][["BP"]]

#moduleGO_res = moduleGO(genes=ddcor_res_GOI$Gene1)

#moduleGO_df = extractModuleGo(moduleGO_res)

### PLOTTING ###

#plotCors(inputMat = vs_counts, design = design_mat, compare = c("GBM", "Control"), geneA = geneOfInterest, geneB = "KIF14")


sig_gain_df_list <- ddcorGO_res[["enrichment_significant_gain_of_correlation_genes"]]
sig_loss_df_list <- ddcorGO_res[["enrichment_significant_loss_of_correlation_genes"]]

labels = c(expression(paste("Ctrl ", rho ," > GBM ", rho)), "GO Term", expression(paste("Ctrl ", rho ," < GBM ", rho)))
#labels = c(expression(paste("Ctrl MER57/MER44B ", rho ," > GBM MER57/MER44B ", rho)), "GO Term Name", expression(paste("Ctrl MER57/MER44B ", rho ," < GBM MER57/MER44B ", rho)))

#plotGOOneGroup(sig_gain_df_list)

#plotGOTwoGroups(sig_gain_df_list,sig_loss_df_list,labels=labels,GOTermTypes =c("BP","CC","MF"))
plotGOTwoGroups(sig_gain_df_list,sig_loss_df_list,labels=labels,GOTermTypes =c("BP")) +
  ggtitle("1q22-correlating GO (Controls)")

# dfList <- sig_gain_df_list
# 
# 
# vscts_top =  filterGenes(vs_counts, 
#                             filterTypes = c("central", "dispersion"), filterCentralPercentile = 0.75, 
#                             filterDispersionPercentile = 0.99)
# 
# ddcor_res = ddcorAll(inputMat = vscts_top, corrType = "spearman", design = design_mat,
#                      compare = c("GBM", "Control"),
#                      adjust = "BH", heatmapPlot = TRUE, nPerm = 0, nPairs = "all")