#########################
#Takes in GTF file specifying transcripts over which to sample
#Creates distribution of exon numbers (for now, assume all exons are equally sized)
#Separately for protein-coding versus noncoding transcripts
#
#
#Takes in another file specifying proviruses that were exonized in annotated genes
#For now, just assume each provirus is the same size and equally likely to be 
#in a given exon
#
#Compares the distribution recovered from randomly picking exons from the
#transcript distribution versus the observed distribution

library(tidyverse)
library(ggplot2)
#source("/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/scripts/SummarySE.R")
#source("/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/scripts/SummarySEwithin.R")
#source("/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/scripts/NormDataWithin.R")
source("C:/Users/worst/iCloudDrive/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/scripts/SummarySE.R")
source("C:/Users/worst/iCloudDrive/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/scripts/SummarySEwithin.R")
source("C:/Users/worst/iCloudDrive/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/scripts/NormDataWithin.R")


options(digits = 20)  


#Base_dir <- "/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis"
Base_dir <- "C:/Users/worst/iCloudDrive/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis"
results_dir = file.path(Base_dir, 'results')

### provirus GTF path ###
#Provirus_gtf_path <- file.path(Base_dir,"Annotations","112922_HERVs_Merged_HiSensitivty_LS_CytoBand_FINAL.tsv")
Provirus_gtf_path <- file.path(Base_dir,"Annotations","083123_hg38_p14_HERVs_HS_LS_CytoBand.tsv")
###For full RefSeq annotation###
#Ref_trans_path <- file.path(Base_dir,"Annotations","T2T_full_ncbiRefSeq_final.txt")

###For full GENCODE annotation###
#Ref_trans_path <- file.path(Base_dir,"Annotations","gencode.v44.chr_patch_hapl_scaff.annotation.gtf")
Ref_trans_path <- file.path(Base_dir,"Annotations","gencode.v44.annotation.gtf")

###For annotation only including transcripts within which exonized HERVs were found
##Read in counts matrix (previously made)
#load(file.path(results_dir,file='telescope','counts_matrices','raw','hg38','Combined','Combined_GBM.Rdata'))
Genome="hg38"
base_dir    = Base_dir
results_dir = file.path(base_dir, 'results')
raw_count_file_dir  = file.path(results_dir,'telescope', 'counts_matrices', 'raw',Genome,'All_Proviruses')

counts = read.csv(paste(raw_count_file_dir,'Glioblastoma_Telescope_output_hg38_p14_RepeatMasker_HS_LS_Reports_LS.csv',sep='/'))
rownames(counts) <- counts$X
counts <- subset(counts,select = -c(X))

#counts <- as.data.frame(countmat)
coldata = read.csv(paste(results_dir,'runTable_GBM.csv',sep='/'))

##Only use GBM data
just_using_runs <- coldata$run[(coldata$disease_state == 'GBM')]
coldata<- subset(coldata,disease_state=='GBM')
counts <- counts[, just_using_runs]
countmat <- as.matrix(counts)
rownames(countmat) <- rownames(counts)

#Removing ALL SHS runs
just_using_runs <- coldata$run[!(coldata$Dataset == 'SHS')]
coldata<- subset(coldata,Dataset!='SHS')

countmat <- countmat[, just_using_runs]


library(GenomicFeatures)
txdb <- makeTxDbFromGFF(Provirus_gtf_path,format="gtf")
# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum them
exonic.gene.sizes <- sum(width(GenomicRanges::reduce(exons.list.per.gene)))
exonic.gene.sizes <- as.data.frame(exonic.gene.sizes)
exonic.gene.sizes$gene <- rownames(exonic.gene.sizes)

#GenesAndLens <- exonic.gene.sizes[rownames(exonic.gene.sizes) %in% genesofinterest, ]
GenesAndLens <- exonic.gene.sizes
names(GenesAndLens) <- c("length","gene")
GenesAndLens <- t(GenesAndLens)
nms <- colnames(GenesAndLens)
GenesAndLens <- as.integer((GenesAndLens[1,]))
names(GenesAndLens) <- nms

# #Normalize to TPM

Counts_to_tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

countmat_sorted <- countmat[order(rownames(countmat)),]

orig_tpm_countmat <-  Counts_to_tpm(subset(countmat_sorted,rownames(countmat_sorted) != ""),GenesAndLens)

#Filter out any provirus with less than 1 TPM median across all samples
# Calculate the median for each row
medians <- apply(orig_tpm_countmat, 1, median)

# Filter the matrix
filtered_tpm_countmat <- orig_tpm_countmat[medians >= 1, ]

###Now, get only proviruses that are exonized ###
## Pull in the names of all exonized proviruses
Exonized_prvs_path <- file.path(Base_dir,"Provirus_Gene_OLs","overlaps_ALLprvs_exons.csv")

Exonized_prvs_base <- read.table(Exonized_prvs_path,sep=',',header=T)
Exonized_prvs_base <- unique(Exonized_prvs_base[c('Provirus','gene_name','transcript_type','is_first_exon','is_last_exon','exon_pos','relative_orientation')])
rownames(Exonized_prvs_base) <- NULL


all_exonic_pvs_in_ProteinCodings <- subset(Exonized_prvs_base, Exonized_prvs_base$transcript_type == "protein_coding")[,"Provirus"]

all_exonic_pvs_in_PCGs <- subset(Exonized_prvs_base, Exonized_prvs_base$transcript_type == "protein_coding")[,"Provirus"]
num_proviruses_PCGs <-length(unique(all_exonic_pvs_in_PCGs))


### Now obtain reference ProteinCodings and PCGs from the GENCODE annotation
Ref_transcripts <- read.table(Ref_trans_path,sep="\t")

#Split the ninth column
max_subfields <- max(str_count(Ref_transcripts$V9, ";")) + 1
Ref_transcripts <- Ref_transcripts %>%
  separate(V9, into = paste0("info_", 1:max_subfields), sep = ";")


##First, extract exon number
Ref_transcripts <- Ref_transcripts %>% 
  mutate(exon_number = str_extract(info_7, " exon_number \\d+"))

#Ref_transcripts <- Ref_transcripts_chk %>%
#  mutate(exon_number = as.numeric(str_extract(exon_number, "\\d+")))

#Subset out the "gene" entries
Ref_transcripts <- subset(Ref_transcripts,V3 != "gene")

Ref_transcripts$info_2 <- as.character(Ref_transcripts$info_2)

#count(Ref_transcripts,"info_2")

#Next, group by transcript and calculate max exon number, and add column
Ref_transcripts <- Ref_transcripts %>%
  group_by(info_2) %>%
  mutate(max_exons = max(exon_number, na.rm = TRUE)) %>%
  ungroup()

#Finally, add exon position columns
Ref_transcripts <- Ref_transcripts %>%
  mutate(is_first_exon = exon_number == 1,
         is_last_exon = exon_number == max_exons,
         is_intermediate_exon = !(is_first_exon | is_last_exon))


ref_ProteinCodings <- subset(Ref_transcripts,info_5==" transcript_type protein_coding")

##Now, subset these so that only ProteinCodings that had proviruses in them are considered
ref_ProteinCodings <- ref_ProteinCodings %>%
  mutate(gene_name = sub(" gene_name ", "", ref_ProteinCodings$info_4))


#ref_ProteinCodings_chk <- subset(ref_ProteinCodings, ref_ProteinCodings$gene_name %in% unique(Exonized_prvs_base$gene_name))
ref_ProteinCoding_exons <- subset(ref_ProteinCodings,V3=="exon")

#ref_PC_trans <- subset(Ref_transcripts,info_5=="	transcript_type protein_coding")
#ref_PC_trans_exons <- subset(ref_PC_trans,info_5=="exon")


ProteinCoding_exon_nums <- ref_ProteinCoding_exons[c("exon_number","max_exons")]
#ProteinCoding_exon_nums_formatted <- cbind(sub("exon_number ","",ProteinCoding_exon_nums$exon_number) , sub("terminal_exon ","", ProteinCoding_exon_nums$max_exons))
#ProteinCoding_exon_nums_formatted <- as.data.frame(cbind(ProteinCoding_exon_nums_formatted[,1], sub(";","", ProteinCoding_exon_nums_formatted[,2])))

names(ProteinCoding_exon_nums) <- c("exon_number","terminal_exon")
ProteinCoding_exon_nums <- as.data.frame(lapply(ProteinCoding_exon_nums, function(x) as.numeric(gsub('exon_number ','', x)) ) )

### Read in actually observed distribution info ###

all_pvs_in_ProteinCodings <- rownames(filtered_tpm_countmat)[rownames(filtered_tpm_countmat) %in% all_exonic_pvs_in_ProteinCodings]
Exonized_prvs_subset <- subset(Exonized_prvs_base, Provirus %in% all_pvs_in_ProteinCodings)

all_pvs_in_ProteinCodings_df <- as.data.frame(cbind(Exonized_prvs_subset$Provirus,Exonized_prvs_subset$exon_pos))

unique_pvs_in_ProteinCodings <- unique(all_pvs_in_ProteinCodings_df)

num_proviruses_in_ProteinCoding <-length(unique_pvs_in_ProteinCodings$V1)

n_dups <- nrow(subset(as.data.frame(table( all_pvs_in_ProteinCodings_df$V2 )),Freq>1))

observed_PRV_dist <- as.data.frame(table(all_pvs_in_ProteinCodings_df$V2,useNA='no'))
observed_PRV_dist <- observed_PRV_dist[,2]
obs_counts_nms <- c('three_prm_terminal_exons',
                    'five_prm_terminal_exons',
                    'single_exon_trans',
                    'nonterminal_exons')[order(c(2,1,4,3))]
names(observed_PRV_dist) <- obs_counts_nms
observed_PRV_dist <- as.data.frame(observed_PRV_dist)

####Begin Monte Carlo####
###Get a random exon,check whether it's the 5' terminal or 3' terminal

###ProteinCodings###
MC_results <- setNames(as.list(c(0,0,0,0)), c('three_prm_terminal_exons',
                                              'five_prm_terminal_exons',
                                              'single_exon_trans',
                                              'nonterminal_exons'))
num_iterations=5000
df_MC_results=data.frame(matrix(ncol = 4, nrow = 0))
colnames(df_MC_results)=names(MC_results)
for(bs in 1:num_iterations){
  
  MC_results <- setNames(as.list(c(0,0,0,0)), c('three_prm_terminal_exons',
                                                'five_prm_terminal_exons',
                                                'single_exon_trans',
                                                'nonterminal_exons'))
  
  for (i in 1:num_proviruses_in_ProteinCoding) {
    random_index=round(runif(1,min=1,max=num_proviruses_in_ProteinCoding))
    selected_exon <- ProteinCoding_exon_nums[random_index,]
    if(selected_exon$terminal_exon == 1){
      MC_results[["single_exon_trans"]]=MC_results[["single_exon_trans"]] + 1
    }
    else{
      if(selected_exon$exon_number == 1)
        MC_results[["five_prm_terminal_exons"]]=MC_results[["five_prm_terminal_exons"]] + 1
      else if(selected_exon$exon_number == selected_exon$terminal_exon)
        MC_results[["three_prm_terminal_exons"]]=MC_results[["three_prm_terminal_exons"]] + 1
      else
        MC_results[["nonterminal_exons"]]=MC_results[["nonterminal_exons"]] + 1
    }
  }
  to_append <- as.data.frame(MC_results)
  df_MC_results <- rbind(df_MC_results, to_append)
  
}

df_MC_results <- as.data.frame(dplyr::mutate_all(as.data.frame(df_MC_results), function(x) as.numeric(x)))



# Convert it to long format
library(reshape2)

df_MC_results$iteration <- as.numeric(rownames(df_MC_results))
data_long <- melt(data=df_MC_results, id.var="iteration",
                  measure.vars=colnames(df_MC_results),
                  variable.name="Exon_Location")
names(data_long)[names(data_long)=="value"] <- "count"



df_MC_results_summary <- summarySEwithin(data_long,measurevar="count",withinvars="Exon_Location",
                                         idvar="iteration",conf.interval=.95)

df_MC_results_summary  <- df_MC_results_summary[-c(2),]
rownames(df_MC_results_summary) <- NULL


#### FIXING SUMMARY ####

df_MC_results_summary[3,"sd"]<-sd(df_MC_results$single_exon_trans)
df_MC_results_summary[1,"sd"]<-sd(df_MC_results$five_prm_terminal_exons)
df_MC_results_summary[2,"sd"]<-sd(df_MC_results$nonterminal_exons)
df_MC_results_summary[4,"sd"]<-sd(df_MC_results$three_prm_terminal_exons)

q_single_exon <- quantile(df_MC_results$single_exon_trans,probs=c(0.025,0.975))
df_MC_results_summary[3,"ci_hi"]<-q_single_exon[2]
df_MC_results_summary[3,"ci_lo"]<-q_single_exon[1]

q_fpte <- quantile(df_MC_results$five_prm_terminal_exons,probs=c(0.025,0.975))
df_MC_results_summary[1,"ci_hi"]<-q_fpte[2]
df_MC_results_summary[1,"ci_lo"]<-q_fpte[1]

q_nte <- quantile(df_MC_results$nonterminal_exons,probs=c(0.025,0.975))
df_MC_results_summary[2,"ci_hi"]<-q_nte[2]
df_MC_results_summary[2,"ci_lo"]<-q_nte[1]

q_tpte <- quantile(df_MC_results$three_prm_terminal_exons,probs=c(0.025,0.975))
df_MC_results_summary[4,"ci_hi"]<-q_tpte[2]
df_MC_results_summary[4,"ci_lo"]<-q_tpte[1]


for (x in 1:4){
  df_MC_results_summary[x,"se"] <- (df_MC_results_summary[x,"sd"])^2/sqrt(df_MC_results_summary[x,"N"])
}


########################

observed_PRV_dist$Exon_Location <- rownames(observed_PRV_dist)
observed_PRV_dist$ObservedOrSimulated <- "Observed"
observed_PRV_dist$count <- observed_PRV_dist$observed_PRV_dist
observed_PRV_dist$counts <- NULL

df_MC_results_summary$ObservedOrSimulated <- "Simulated"

full_df <- rbind.fill(df_MC_results_summary,observed_PRV_dist)

##To plot the observed data
# ggplot(as.data.frame(observed_PRV_dist), aes(x=factor(rownames(observed_PRV_dist),level=c('three_prm_terminal_exons',
#                                                                                           'five_prm_terminal_exons',
#                                                                                           'single_exon_trans',
#                                                                                           'nonterminal_exons')), y=counts)) + 
#   geom_bar(position=position_dodge(), stat="identity",
#            colour="black", # Use black outlines,
#            size=.3) +      # Thinner lines
#   xlab("Exon Location") +
#   ylab("Count") +
#   ggtitle("Observed HERV distribution in ProteinCodings") +
#   theme_bw()

##To plot the observed AND the simulated ##


full_df$ObservedOrSimulated <- factor(full_df$ObservedOrSimulated, levels = c("Simulated", "Observed"))

ggplot(as.data.frame(full_df), aes(x=Exon_Location, y=count, fill=ObservedOrSimulated)) + 
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  geom_errorbar(aes(ymin=ci_lo, ymax=ci_hi),
                size=1.0,    # Thick lines
                width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(name="Observed\nor\nSimulated",
                    values = c("red", "blue"),
                    limits = c("Simulated","Observed")) +
  xlab("Exon Location") +
  ylab("Count") +
  ggtitle("Simulated HERV distribution in mRNAs") +
  theme_bw() +
  theme(legend.title = element_text(size = 14, 
                                    face = "bold", 
                                    hjust = 0.5, 
                                    vjust = 1.2, 
                                    lineheight = 0.8),
        legend.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 16, color = "black"), # Larger and darker axis text
        plot.title = element_text(size = 28, hjust = 0.5)) +
  labs(fill = "") +
  scale_x_discrete(labels = c("5' Terminal", "Nonterminal", "Single", "3' Terminal")) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) # Ensure each line has a label

### GET p VALUES ###

df_MC_results$three_prm_terminal_exons_withinCI <- NA
df_MC_results$five_prm_terminal_exons_withinCI <- NA
df_MC_results$single_exon_withinCI <- NA
df_MC_results$nonterminal_exons_withinCI <- NA

for(i in 1:nrow(df_MC_results)){
  
  ##3 prime terminal exons (seem to be enriched in observed data)
  if(observed_PRV_dist[1,] <= df_MC_results[i,"three_prm_terminal_exons"])
  {
    df_MC_results[i,6]=1
  }
  else
    df_MC_results[i,6]=0
  
  ##5 prime terminal exons (seem actually to be depleted in observed data)
  if(observed_PRV_dist[3,] >= df_MC_results[i,"five_prm_terminal_exons"])
  {
    df_MC_results[i,7]=1
  }
  else
    df_MC_results[i,7]=0
  
  ##single exons (seem to be enriched in observed data)
  if(observed_PRV_dist[2,] <= df_MC_results[i,"single_exon_trans"])
  {
    df_MC_results[i,8]=1
  }
  else
    df_MC_results[i,8]=0
  
  ## nonterminal exons (seem to be depleted in observed data)
  if(observed_PRV_dist[4,] >= df_MC_results[i,"nonterminal_exons"])
  {
    df_MC_results[i,9]=1
  }
  else
    df_MC_results[i,9]=0
}


pval_3pmtermexons <- (sum(df_MC_results$three_prm_terminal_exons_withinCI)+1)/(nrow(df_MC_results)+1)
pval_3pmtermexons

pval_5pmtermexons <- (sum(df_MC_results$five_prm_terminal_exons_withinCI)+1)/(nrow(df_MC_results)+1)
pval_5pmtermexons

pval_single_exons <- (sum(df_MC_results$single_exon_withinCI)+1)/(nrow(df_MC_results)+1)
pval_single_exons

pval_nonterm_exons <- (sum(df_MC_results$nonterminal_exons_withinCI)+1)/(nrow(df_MC_results)+1)
pval_nonterm_exons

options(digits = 5)  


### FOR PC GENES ###
### Read in actually observed distribution info ###
overexp_pvs_unique_pcs_path <- file.path(Exonized_prvs_path,'Overexpressed',
                                         'overexp_pvs_PCgenes.txt')
overexp_pvs_unique_pcs <- read.table(overexp_pvs_unique_pcs_path,sep=" ",fill=TRUE)

overexp_pvs_unique_pcs$UTR_Location <- paste(overexp_pvs_unique_pcs$V2,overexp_pvs_unique_pcs$V3,overexp_pvs_unique_pcs$V4)
overexp_pvs_unique_pcs <- overexp_pvs_unique_pcs[,c(1,5)]

underexp_pvs_unique_pcs_path <- file.path(Exonized_prvs_path,'Underexpressed',
                                          'underexp_pvs_PCgenes.txt')
underexp_pvs_unique_pcs <- read.table(underexp_pvs_unique_pcs_path,sep=" ",fill=TRUE)

underexp_pvs_unique_pcs$UTR_Location <- paste(underexp_pvs_unique_pcs$V2,underexp_pvs_unique_pcs$V3,underexp_pvs_unique_pcs$V4,underexp_pvs_unique_pcs$V5,underexp_pvs_unique_pcs$V6)
underexp_pvs_unique_pcs <- underexp_pvs_unique_pcs[,c(1,7)]

all_pvs_unique_pcs <- rbind(underexp_pvs_unique_pcs,overexp_pvs_unique_pcs)
n_dups <- nrow(subset(as.data.frame(table(all_pvs_unique_pcs$V1)),Freq>1))

cts <- table(all_pvs_unique_pcs$UTR_Location,useNA='no')

num_3_prime_UTR <- cts[3] + cts[4]
num_5_prime_UTR <- cts[5] + cts[6]
num_both_UTRs <- cts[7]

#obs_counts_nms <- c('three_prm_terminal_exons',
#                    'five_prm_and_three_prm',
#                    'three_prm_terminal_exons',
#                    'CDS')[order(c(1,3,2,4))]
#observed_PRV_dist <- data.frame(as.list(c(num_3_prime_UTR,num_5_prime_UTR,num_both_UTRs)))
#rownames(observed_PRV_dist) <- obs_counts_nms




###BOX AND WHISKER###

# Define the desired exon positions in order
exon_order <- c("five_prm_terminal_exons", "nonterminal_exons", "single_exon_trans", "three_prm_terminal_exons")

# Filter data_long to only include the specified exon positions
filtered_data_long <- data_long %>%
  filter(Exon_Location %in% exon_order) %>%
  mutate(Exon_Location = factor(Exon_Location, levels = exon_order))

# Ensure observed_PRV_dist follows the same exon ordering
observed_PRV_dist <- observed_PRV_dist %>%
  filter(Exon_Location %in% exon_order) %>%
  mutate(Exon_Location = factor(Exon_Location, levels = exon_order))

# Create the plot
ggplot(filtered_data_long, aes(x = Exon_Location, y = count)) +
  geom_boxplot(aes(fill = "Simulated"), alpha = 0.7, outlier.shape = NA) +  # Boxplot for simulated data
  geom_segment(data = observed_PRV_dist, 
               aes(x = as.numeric(Exon_Location)-0.2, xend = as.numeric(Exon_Location)+0.2,
                   y = count, yend = count, color = "Observed"), 
               size = 1.8) +  # Ensure observed lines are thick and visible
  scale_x_discrete(labels = c("5' Terminal", "Nonterminal", "Single", "3' Terminal")) +  
  guides(x.sec = "axis", y.sec = "axis") +
  scale_y_break(c(500, 2300), space = 2) +  # Break in y-axis
  scale_fill_manual(name = "Legend", values = c("Simulated" = "#E63946")) +  # Vivid red for simulated data
  scale_color_manual(name = "Legend", values = c("Observed" = "#1D3557")) +  # Blue for observed median lines
  labs(title = NULL,
       x = "Exon Location",
       y = "Count") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),  # Remove "fill" and "colour" from the legend
    legend.text = element_text(size = 14),  # Increase legend text size
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.title = element_text(size = 16),  # Increase axis title size
    axis.text = element_text(size = 14, color = "black"),  # Increase tick labels size and darken color
    axis.ticks = element_line(size = 1.2),  # Make axis ticks more visible
    axis.line = element_line(size = 1.5, color = "black"),  # Thick x and y axis lines
    panel.grid.major = element_line(color = "gray80"),  # Improve grid contrast
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x.top = element_blank(),  # Remove duplicated x-axis labels at the top
    axis.ticks.x.top = element_blank(),  # Remove top x-axis ticks
    panel.grid.major.y = element_blank()  # Remove y-grid to leave space for break symbol
  )
