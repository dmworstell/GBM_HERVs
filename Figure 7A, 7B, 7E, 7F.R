###############
#Takes in combined counts data from telescope
#Takes in a file with metadata including variables to compare
#Analyzes 1q22 survival including taking average of counts when there are
#Multiple patients

#Daniel Murimi-Worstell
#6/4/24

#Edited:

###############

#Load packages
suppressPackageStartupMessages(library(genefilter))
suppressPackageStartupMessages(library(apeglm))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(PCAtools))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(EBSeq))


#Load in data
base_dir    = '/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis'
results_dir = file.path(base_dir, 'results')
Genome="hg38"
raw_count_file_dir  = file.path(results_dir,'telescope', 'counts_matrices', 'raw',Genome,'HML2_Genes')
counts = read.csv(paste(raw_count_file_dir,'Glioblastoma_Telescope_output_HML2_GENES.csv',sep='/'))

###Get the HERV sheet to combine###
raw_count_file_dir_HERV  = file.path(results_dir,'telescope', 'counts_matrices', 'raw',Genome,'All_Proviruses')
GBM_counts_HERV = read.csv(paste(raw_count_file_dir_HERV,'Glioblastoma_Telescope_output_hg38_p14_RepeatMasker_HS_LS_Reports_LS.csv',sep='/'))
counts <- rbind(counts,GBM_counts_HERV)

rm(GBM_counts_HERV)

#Convert counts to matrix
countmat = round(data.matrix(counts[,-1]))
rownames(countmat) = counts$X
rm(counts)

#######

#Read in column data
coldata = read.csv(paste(results_dir,'runTable_GBM.csv',sep='/'))

#Remove rows that weren't run through telescope
coldata <- coldata[coldata$run %in% colnames(countmat), ] 
rownames(coldata) <- NULL

#Removing SHS runs
just_using_runs <- coldata$run[!(coldata$Dataset == 'SHS')]
coldata<- subset(coldata,Dataset!='SHS')
countmat <- countmat[, just_using_runs]

#Go to the correct directory
setwd('/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/figures')

#Filter rows with median counts under 2
keep = which(rowMeans(countmat) > 2)
countmat <- countmat[keep,]

####PREPROCESSING####
#If primary and recurrent samples from the same subject, remove recurrent

# Step 1: Identify non-unique patients
non_unique_patients <- coldata %>%
  group_by(patient) %>%
  filter(n() > 1)

# Step 2: For non-unique patients, keep rows where recurrence is "Primary"
primary_rows <- non_unique_patients %>%
  filter(recurrence == "Primary")

# Step 3: Combine the primary rows with the rest of the dataframe
# Remove the non-unique patients from the original dataframe
unique_patients <- coldata %>%
  filter(!patient %in% non_unique_patients$patient)

# Add the primary rows back to the dataframe
coldata_filtered <- bind_rows(unique_patients, primary_rows)

# Step 4: If there are still multiple rows for any patients, keep only the first one
coldata_final <- coldata_filtered %>%
  group_by(patient) %>%
  slice(1) %>%
  ungroup()

#Extract the run values from coldata_final
selected_runs <- coldata_final$run

coldata <- coldata_final
rm(coldata_final,coldata_filtered)

# Subset the columns of countmat to include only the selected runs
countmat_filtered <- countmat[, selected_runs, drop = FALSE]

countmat <- countmat_filtered
rm(countmat_filtered)
#Otherwise, if more than one primary, take the first one


#Generate Normalized counts matrix
Norm_counts_mat <- GetNormalizedMat(countmat, MedianNorm(countmat))
rm(countmat)
scaledata <- t(scale(t(Norm_counts_mat)))
rm(Norm_counts_mat)



## BY SPECIFIC PROVIRUS
  
## ERVK-7, HML-2 1q22, HML-K102
#prv_choice <- "1q22_chr1:155627723-155634872_HERVK-int"
#prv_choice <- "HML-2_1q22"

#prv_label <- "HML-2 1q22"
#prv_label <- "LTR12F"
#prv_label <- "2q12.1 MER57-int"
#prv_label <- "CARKD"

## MER57-int provirus
#prv_choice <- "2q12.1_chr2:103477983-103481010_MER57-int"

##ERVK3-1
prv_choice <- "19q13.43_chr19:58305671-58315267_LTR3"
prv_label <- "ERVK-3"

##ID'd in screen##
#prv_choice <- "13q34_chr13:110613521-110614348_LTR12F"

#Associated gene#
#prv_choice <- "CARKD"

sp_counts_unaveraged <- as.data.frame(scaledata[prv_choice,])

survival_unaveraged <- select(coldata,c("patient","OS","PFS","Vital.Status.1.Dead"))
survival_unaveraged_unique_rows <- distinct(survival_unaveraged)

survival_unaveraged_with_runs <- select(coldata,c("run","patient","OS","PFS","Vital.Status.1.Dead"))

##Variance stabilized, transformed counts
#sp_counts_unaveraged <- as.data.frame(vs_counts[prv_choice,])

colnames(sp_counts_unaveraged) <- "norm_counts"
sp_counts_with_survival <- cbind(sp_counts_unaveraged,survival_unaveraged)


##DON'T AVERAGE, take primary sample 
median_counts_per_pt <- aggregate(norm_counts~patient,data=sp_counts_with_survival, FUN=mean)


median_counts_per_pt <- aggregate(norm_counts~patient,data=sp_counts_with_survival, FUN=median)
median_counts_with_survival_per_pt <- cbind(survival_unaveraged_unique_rows,median_counts_per_pt)

rownames(median_counts_with_survival_per_pt) <- median_counts_with_survival_per_pt$patient

kmc_df <- data.frame(select(median_counts_with_survival_per_pt,c("norm_counts","OS","PFS","Vital.Status.1.Dead")))

kmc_df_noNA <- subset(kmc_df,OS!="UNKNOWN")
kmc_df_noNA$OS <- as.numeric(kmc_df_noNA$OS)
kmc_df_noNA$PFS <- as.numeric(kmc_df_noNA$PFS)
kmc_df_noNA$Vital.Status.1.Dead <- as.numeric(kmc_df_noNA$Vital.Status.1.Dead)
kmc_df_sorted <- kmc_df_noNA[order(kmc_df_noNA$norm_counts),]

os_fit <- survfit(Surv(OS, Vital.Status.1.Dead) ~ 1, data = kmc_df_sorted)
pfs_fit <- survfit(Surv(PFS, Vital.Status.1.Dead) ~ 1, data = kmc_df_sorted)

surv_object <- Surv(time = kmc_df_sorted$OS, event = kmc_df_sorted$Vital.Status.1.Dead)

# create a Kaplan-Meier plot stratified by transformed count quartiles
#kmc_df_sorted$quartile <- cut(kmc_df_sorted$norm_counts, breaks = quantile(kmc_df_sorted$norm_counts, probs = seq(0, 1, 0.25)), labels = FALSE,include.lowest=T)
#legendlabs <- c(paste0(prv_label," Q1"), paste0(prv_label," Q2"), paste0(prv_label," Q3"), paste0(prv_label," Q4"))

# create a Kaplan-Meier plot stratified by transformed count high or low
kmc_df_sorted$quartile <- cut(kmc_df_sorted$norm_counts, breaks = quantile(kmc_df_sorted$norm_counts, probs = c(0,0.5,1)), labels = FALSE,include.lowest=T)
legendlabs <- c(paste0(prv_label," LOW"), paste0(prv_label," HIGH"))

km_fit <- survfit(Surv(OS, Vital.Status.1.Dead) ~ quartile, data = kmc_df_sorted)

ggsurvplot(
  km_fit, 
  data = kmc_df_sorted, 
  pval = TRUE, 
  risk.table = FALSE,  # Remove the number at risk table
  conf.int = TRUE, 
  legend.labs = legendlabs,
  legend.title = element_blank(),  # Remove the strata label in the legend
  ggtheme = theme_minimal(),  # Apply a minimal theme for customization
  font.main = c(20, "bold"),  # Main title font size and style
  font.x = c(18, "bold"),     # X-axis label font size and style
  font.y = c(18, "bold"),     # Y-axis label font size and style
  font.tickslab = c(16, "plain"),  # Axis ticks font size and style
  font.legend = c(13, "plain")  # Legend font size and style
)
####