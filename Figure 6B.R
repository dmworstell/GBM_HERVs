library(tidyverse)
library(DGCA)
library(psych)
library(cluster)
library(reshape2)
suppressPackageStartupMessages(library(apeglm))
suppressPackageStartupMessages(library(DESeq2))
library(ggplot2, quietly = TRUE)
library(GOstats, quietly = TRUE)
library(HGNChelper, quietly = TRUE)
library(org.Hs.eg.db, quietly = TRUE)
library(plotrix, quietly = TRUE)
library(rstatix)

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

herv_of_interest <- "2q12.1_chr2:103477983-103481010_MER57-int"
GBM_counts_HERV_chk <- subset(GBM_counts_HERV,X==herv_of_interest)
counts <- rbind(counts,GBM_counts_HERV_chk)
#counts <- rbind(counts,GBM_counts_HERV)

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

#Just TCGA runs
just_using_runs <- coldata$run[(coldata$study_name == 'TCGA')]
coldata<- subset(coldata,study_name=='TCGA')
countmat <- countmat[, just_using_runs]

### Analyze correlations

library(EBSeq)

Norm_counts_mat_unfiltered <- GetNormalizedMat(countmat, MedianNorm(countmat))

## NON-FILTERED
# Step 1: Calculate the mean for each gene across all samples
prv_means <- rowMeans(Norm_counts_mat_unfiltered)

# Step 2: Calculate quartiles from the vector of means
Q1 <- quantile(prv_means, 0.25)  # 25th percentile, first quartile
Q2 <- quantile(prv_means, 0.50)  # 50th percentile, second quartile (median)
Q3 <- quantile(prv_means, 0.75)  # 75th percentile, third quartile

# Printing the quartiles
print(paste("First Quartile (Q1):", Q1))
print(paste("Second Quartile (Q2, Median):", Q2))
print(paste("Third Quartile (Q3):", Q3))


## FILTERED
keep = which(rowMeans(Norm_counts_mat_unfiltered) > 50)
Norm_counts_mat <- Norm_counts_mat_unfiltered[keep,]

# Get mean expression quartiles

# Step 1: Calculate the mean for each gene across all samples
prv_means <- rowMeans(Norm_counts_mat)

# Step 2: Calculate quartiles from the vector of means
Q1 <- quantile(prv_means, 0.25)  # 25th percentile, first quartile
Q2 <- quantile(prv_means, 0.50)  # 50th percentile, second quartile (median)
Q3 <- quantile(prv_means, 0.75)  # 75th percentile, third quartile

# Printing the quartiles
print(paste("First Quartile Normalized (Q1):", Q1))
print(paste("Second Quartile Normalized (Q2, Median):", Q2))
print(paste("Third Quartile Normalized (Q3):", Q3))


##Assess normality##
# Extract expression data for 2q12.1_chr2:103477983-103481010_MER57-int
prv_expression <- Norm_counts_mat["2q12.1_chr2:103477983-103481010_MER57-int", ]

# Generate a Q-Q plot
qqnorm(prv_expression)
qqline(prv_expression, col = "red")  # Adds a reference line

# Perform Shapiro-Wilk test
shapiro_test_result <- shapiro.test(prv_expression)
print(shapiro_test_result)


# Obtain the central region corresponding to Normality
MER57int_expression <- Norm_counts_mat["2q12.1_chr2:103477983-103481010_MER57-int", ]

lower_quantile_value <- quantile(MER57int_expression, probs = 0.25) # Approximation for theoretical quantile -1
upper_quantile_value <- quantile(MER57int_expression, probs = 0.75) # Approximation for theoretical quantile 1

# Identify samples within these quantile values
samples_within_range <- which(MER57int_expression >= lower_quantile_value & MER57int_expression <= upper_quantile_value)

# Subset Norm_counts_mat to retain only those samples
subset_Norm_counts_mat <- Norm_counts_mat[, samples_within_range]

#Assess normality
prv_expression <- subset_Norm_counts_mat["2q12.1_chr2:103477983-103481010_MER57-int", ]
# Generate a Q-Q plot
qqnorm(prv_expression)
qqline(prv_expression, col = "red")  # Adds a reference line

# Perform Shapiro-Wilk test
shapiro_test_result <- shapiro.test(prv_expression)
print(shapiro_test_result)

countmat_filtered_cd <- Norm_counts_mat

# identifying HERVs based on underscores anywhere in the name or HML2
is_herv <- grepl("HML-2|HML2.LTR", rownames(countmat_filtered_cd)) | sapply(strsplit(rownames(countmat_filtered_cd), "_"), function(x) length(x) - 1) >= 2

# Keep only HERVs of interest and all genes
filtered_countmat <- countmat_filtered_cd[rownames(countmat_filtered_cd) %in% herv_of_interest | !is_herv, ]

## Generate correlation matrices ##

# Extract only gene rows from the count matrix
genes_only <- filtered_countmat[!(rownames(filtered_countmat) %in% herv_of_interest),]

variance_filter <- apply(genes_only, 1, var) != 0
genes_only <- genes_only[variance_filter, ]

# Initialize a list to store significantly correlated genes for each HERV
significant_correlations <- list()

  herv <- herv_of_interest
  # Extract the expression data for the current HERV
  herv_expression <- as.numeric(filtered_countmat[rownames(filtered_countmat) == herv, ])
  
  # Initialize vectors to store p-values and gene names
  p_values <- numeric(nrow(genes_only))
  correlation_coefficients <- numeric(nrow(genes_only))
  gene_names <- rownames(genes_only)
  
  # Perform correlation analysis for each gene
  for (i in 1:nrow(genes_only)) {
    gene_expression <- as.numeric(genes_only[i, ])
    
    # Use cor.test to get the p-value
    test_result <- cor.test(herv_expression, gene_expression, method = "spearman", use = "pairwise.complete.obs")
    
    # Store the p-value
    p_values[i] <- test_result$p.value
    correlation_coefficients[i] <- test_result$estimate
  }
  
  # Apply the Benjamini-Hochberg procedure to adjust p-values for multiple testing
  p_adjusted <- p.adjust(p_values, method = "BH")
  
  # Select genes with adjusted p-value < 0.05
  significant_indices <- which(p_adjusted < 0.05)
  significant_genes <- gene_names[p_adjusted < 0.05]
  
  # Create a data frame with gene names, correlation coefficients, and p-values
  significant_data <- data.frame(
    gene = significant_genes,
    correlation_coefficient = correlation_coefficients[significant_indices],
    p_value = p_values[significant_indices],
    p_adjusted = p_adjusted[significant_indices]
  )
  
  # Store the results for the current HERV
  significant_correlations[[herv]] <- significant_genes

library(clusterProfiler)

  # Define a custom theme
  custom_theme <- function(font.size = 12, label.font.size = 10) {
    theme_classic(base_size = font.size) +
      theme(
        axis.title.x = element_text(size = font.size),
        axis.title.y = element_text(size = font.size),
        axis.text.x = element_text(size = font.size),
        axis.text.y = element_text(size = label.font.size),  # Adjust label font size
        legend.title = element_text(size = font.size),
        legend.text = element_text(size = font.size),
        plot.title = element_text(size = font.size, hjust = 0.5),
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "right",
        panel.grid.major = element_line(color = "grey80", size = 0.5),
        panel.grid.minor = element_line(color = "grey90", size = 0.25),
        legend.box.margin = margin(0, 0, 0, 12)  # Adjust margin to separate legend from plot
      )
  }
  
  # Define a function to format labels
  .label_format <- function(label_format = 30) {
    function(x) {
      sapply(x, function(y) {
        # Capitalize the first letter of each word and wrap the text
        wrapped_label <- stringr::str_wrap(tools::toTitleCase(y), width = label_format)
        return(wrapped_label)
      })
    }
  }
  
  # Define the custom dotplot function
  custom_dotplot_enrichResult <- function(object, x = "geneRatio", color = "p.adjust",
                                          showCategory=10, size=NULL, split = NULL,
                                          font.size=12, title = "", orderBy="x",
                                          label_format = 30, decreasing=TRUE) {
    if (x == "geneRatio" || x == "GeneRatio") {
      x <- "GeneRatio"
      if (is.null(size))
        size <- "Count"
    } else if (x == "count" || x == "Count") {
      x <- "Count"
      if (is.null(size))
        size <- "GeneRatio"
    } else if (is(x, "formula")) {
      x <- as.character(x)[2]
      if (is.null(size))
        size <- "Count"
    } else {
      if (is.null(size))
        size  <- "Count"
    }
    
    if (inherits(object, c("enrichResultList", "gseaResultList"))) {
      ldf <- lapply(object, fortify, showCategory=showCategory, split=split)
      df <- dplyr::bind_rows(ldf, .id="category")
      df$category <- factor(df$category, levels=names(object))
    } else {
      df <- fortify(object, showCategory = showCategory, split=split)
    }
    
    if (orderBy !=  'x' && !orderBy %in% colnames(df)) {
      message('wrong orderBy parameter; set to default `orderBy = "x"`')
      orderBy <- "x"
    }
    
    if (orderBy == "x") {
      df <- dplyr::mutate(df, x = eval(parse(text=x)))
    }
    
    label_func <- .label_format(label_format)
    
    idx <- order(df[[orderBy]], decreasing = decreasing)
    
    df$Description <- factor(df$Description,
                             levels=rev(unique(df$Description[idx])))
    
    p <- ggplot(df, aes(x=.data[[x]], y=.data[["Description"]], 
                        size=.data[[size]])) +
      geom_point(color = "black", fill = "red", shape = 21, stroke = 0.5) +  # Black border, red fill
      scale_y_discrete(labels = label_func) +
      ylab(NULL) + ggtitle(title) + custom_theme(font.size) +
      guides(size = guide_legend(override.aes = list(shape = 21, fill = "red", color = "black")))  # Legend styling
    
    if (size == "Count") {
      size_break <- pretty(df[[size]], n=4)
      p <- p + scale_size(range=c(3, 8), breaks = size_break)
    } else {
      p <- p + scale_size(range=c(3, 8))
    }
    
    class(p) <- c("enrichplotDot", class(p))
    return(p)
  }
  
  
  
  library(clusterProfiler)
  
  # Define the function to perform GO analysis
  perform_go_analysis <- function(gene_list, orgDb, showCategory = 5, output_path) {
    # Ensure the gene list is a character vector
    gene_list <- as.character(gene_list)
    
    # Prepare gene list for enrichment analysis
    ego <- enrichGO(gene = gene_list,
                    OrgDb = orgDb,
                    keyType = "SYMBOL",
                    ont = "BP",  # CC = Cellular Component, BP = Biological Process, MF = Molecular Function
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)
    
    # Plot GO enrichment result
    if (length(ego) > 0) {  # Check if there are any enriched GO terms
      # Use the custom dotplot function
      plot_enrichment <- custom_dotplot_enrichResult(ego, showCategory = showCategory, label_format = 20) + 
        xlab("Gene Ratio") +
        custom_theme(font.size = 12, label.font.size = 12) +  # Use the custom theme with larger text
        scale_size(name = "Gene Counts", range = c(0.5, 10), breaks = c(150, 200, 250, 350), limits = c(100, 350)) +
        ggtitle("MER56-int GO (GBM)")
      
      # Print the plot
      print(plot_enrichment)
      
      # Extract genes from the plotted GO terms
      plotted_terms <- plot_enrichment$data$ID
      # Create a nested list for genes in each plotted GO term
      genes_in_plotted_terms <- lapply(plotted_terms, function(term_id) {
        # Get the gene list for the current GO term
        gene_list <- ego@result[ego@result$ID == term_id, "geneID"]
        # Split the gene list by '/'
        genes <- unlist(strsplit(gene_list, "/"))
        return(genes)
      })
      # Name the list elements with the GO term IDs
      names(genes_in_plotted_terms) <- plotted_terms
      
      # Initialize a list to store data frames for each GO term
      go_term_data_frames <- list()
      
      # Create a data frame for each GO term
      for (term_id in plotted_terms) {
        genes <- genes_in_plotted_terms[[term_id]]
        
        # Find the intersection with significant genes and get correlation data
        significant_genes_in_term <- significant_data[significant_data$gene %in% genes, ]
        
        # Create a data frame with genes, correlation coefficients, and p-values
        term_df <- data.frame(
          gene = significant_genes_in_term$gene,
          correlation_coefficient = significant_genes_in_term$correlation_coefficient,
          p_value = significant_genes_in_term$p_value
        )
        # Store the data frame in the list
        go_term_data_frames[[term_id]] <- term_df
      }
      
      # Find the maximum number of genes in any GO term
      max_rows <- max(sapply(go_term_data_frames, nrow))
      
      # Function to pad a data frame with NA rows to match the maximum number of rows
      pad_with_na <- function(df, max_rows) {
        if (nrow(df) < max_rows) {
          # Calculate the number of rows to add
          rows_to_add <- max_rows - nrow(df)
          # Create a data frame with NA values
          na_df <- data.frame(matrix(NA, nrow = rows_to_add, ncol = ncol(df)))
          colnames(na_df) <- colnames(df)
          # Bind the NA rows to the original data frame
          df <- bind_rows(df, na_df)
        }
        return(df)
      }
      
      go_term_descriptions <- ego@result[ego@result$ID %in% plotted_terms, "Description"]
      
      # Pad each data frame to have the same number of rows
      padded_data_frames <- lapply(go_term_data_frames, pad_with_na, max_rows = max_rows)
      
      # Combine all data frames into one, with each GO term having its own set of columns
      combined_df <- do.call(cbind, lapply(names(padded_data_frames), function(term_id) {
        term_df <- padded_data_frames[[term_id]]
        # Rename columns to include the GO term ID
        colnames(term_df) <- paste(term_id, colnames(term_df), sep = "_")
        return(term_df)
      }))
      
      # Add the GO term descriptions as the first row
      description_row <- unlist(lapply(names(go_term_data_frames), function(term_id) {
        description <- go_term_descriptions[plotted_terms == term_id]
        # Repeat the description for each column of the term
        rep(description, 3)  # Assuming 3 columns per GO term: gene, correlation, p-value
      }))
      # Create a data frame for the description row
      description_df <- as.data.frame(t(description_row), stringsAsFactors = FALSE)
      colnames(description_df) <- colnames(combined_df)
      
      # Combine the description row with the combined data frame
      final_df <- rbind(description_df, combined_df)
      
      # Write the final data frame to a CSV file
      # write.csv(final_df, file = output_path, row.names = FALSE)
      
      # Return the final data frame
      return(final_df)
      
    } else {
      cat("No significant GO terms found for this set.\n")
    }
  }
  
  # Specify the output path for the CSV file
  output_path <- "/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Manuscripts/2024/Supplemental/Supplemental Data Sheets/final_df.csv"
  
  # Perform GO analysis with a specified number of categories to show and save the result
  final_df <- perform_go_analysis(significant_correlations[[1]], org.Hs.eg.db, showCategory = 5, output_path = output_path)
  # Specify the output path for the CSV file
  output_path <- "/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Manuscripts/2024/Supplemental/Supplemental Data Sheets/final_df.csv"
  
  # Perform GO analysis with a specified number of categories to show and save the result
  final_df <- perform_go_analysis(significant_correlations[[1]], org.Hs.eg.db, showCategory = 5, output_path = output_path)