library(VennDiagram)
library("gridExtra")                # Load gridExtra package
library(ggplot2)

top_prvs_base_dir <- "/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing"

bulk_base_dir <- file.path(top_prvs_base_dir,"TCGA_Analysis/top_proviruses/hg38/Telescope/Combined/TCGA_Alone")
sc_base_dir <- file.path(top_prvs_base_dir,"scRNAseq_Analysis/results/top_proviruses/Matched_With_Bulk/Individual_Proviruses")


#Read in bulk telescope unfiltered results for prvs, hg38
proviruses_bulk <- rownames(read.csv("/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/top_proviruses/090623_results_GBM_vs_Controls_top_proviruses_withSense.csv",row.names=1))

### Venn diagram of Bulk RNA-seq versus scRNA-seq proviruses
#(Neoplastic vs Astrocytes, OPCs or Oligodendrocytes, and Neurons vs these cell types)
sc_ummatched_base_dir <- file.path(top_prvs_base_dir,"scRNAseq_Analysis/results/top_proviruses/Contrast_Method")

um_proviruses_sc_Neoplastic_Astrocyte <- rownames(read.csv(file.path(sc_ummatched_base_dir,
                                                                  "top_proviruses_Neoplastic_vs_Astrocyte_ContrastMethod.csv"),row.names=1))
um_proviruses_sc_Neoplastic_OPC <- rownames(read.csv(file.path(sc_ummatched_base_dir,
                                                            "top_proviruses_Neoplastic_vs_OPC_ContrastMethod.csv"),row.names=1))
um_proviruses_sc_Neoplastic_Oligodendrocyte <- rownames(read.csv(file.path(sc_ummatched_base_dir,
                                                                        "top_proviruses_Neoplastic_vs_Oligodendrocyte_ContrastMethod.csv"),row.names=1))

um_proviruses_sc_Astrocyte_Neuron <- rownames(read.csv(file.path(sc_ummatched_base_dir,
                                                              "top_proviruses_Astrocyte_vs_Neuron_ContrastMethod.csv"),row.names=1))
um_proviruses_sc_OPC_Neuron <- rownames(read.csv(file.path(sc_ummatched_base_dir,
                                                        "top_proviruses_OPC_vs_Neuron_ContrastMethod.csv"),row.names=1))
um_proviruses_sc_Oligodendrocyte_Neuron <- rownames(read.csv(file.path(sc_ummatched_base_dir,
                                                                    "top_proviruses_Oligodendrocyte_vs_Neuron_ContrastMethod.csv"),row.names=1))

proviruses_sc_Neoplastic_Glial <- unique(c(um_proviruses_sc_Neoplastic_Astrocyte,um_proviruses_sc_Neoplastic_OPC,um_proviruses_sc_Neoplastic_Oligodendrocyte))
proviruses_sc_Glial_Neuron <- unique(c(um_proviruses_sc_Astrocyte_Neuron,um_proviruses_sc_OPC_Neuron,um_proviruses_sc_Oligodendrocyte_Neuron))

proviruses_bulk_AND_proviruses_sc_Neoplastic_Glial <- intersect(proviruses_bulk,proviruses_sc_Neoplastic_Glial)
proviruses_sc_Neoplastic_Glial_AND_proviruses_sc_Glial_Neuron <- intersect(proviruses_sc_Neoplastic_Glial,proviruses_sc_Glial_Neuron)
proviruses_bulk_AND_proviruses_sc_Glial_Neuron <- intersect(proviruses_bulk,proviruses_sc_Glial_Neuron)

proviruses_bulk_AND_proviruses_sc_Neoplastic_Glial_AND_proviruses_sc_Glial_Neuron <- intersect(proviruses_bulk_AND_proviruses_sc_Neoplastic_Glial,proviruses_sc_Glial_Neuron)

grid.newpage()
combined.venn.plot <- draw.triple.venn(
  area1 = length(proviruses_bulk),
  area2 = length(proviruses_sc_Neoplastic_Glial),
  area3 = length(proviruses_sc_Glial_Neuron),
  n12 = length(proviruses_bulk_AND_proviruses_sc_Neoplastic_Glial),
  n23 = length(proviruses_sc_Neoplastic_Glial_AND_proviruses_sc_Glial_Neuron),
  n13 = length(proviruses_bulk_AND_proviruses_sc_Glial_Neuron),
  n123 = length(proviruses_bulk_AND_proviruses_sc_Neoplastic_Glial_AND_proviruses_sc_Glial_Neuron),
  #category = c("GBM vs. Controls (Bulk)", "Neoplastic vs. Glial Cells (sc)","Glial Cells vs. Neurons (sc)"),
  fill = c("skyblue", "pink", "green"),
  lty = "blank",
  euler.d=TRUE,
  scaled = TRUE
)

#check which vectors contain specific proviruses
#pv_of_interest <- "1q22_chr1:155626666-155627633_LTR5_Hs"
#pv_of_interest <- "2q12.1_chr2:103477983-103481010_MER57-int"
#for (v in ls()){
#  #print(paste("CHECKING",v))
#  vector <- get(v)
#  if(pv_of_interest %in% vector) {
#    print(v)
#  }
#}