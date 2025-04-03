library(Biostrings)
library(ggseqlogo)
require(ggplot2)

# Define the directories
control_dir <- "/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs/Control_NON_terminal_HERV_fastas/polyA_Contexts/specified_clades"
terminal_dir <- "/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs/terminal_HERV_fastas/polyA_Contexts/specified_clades"

read_fasta_files <- function(clade) {
  # Create file patterns based on the clade
  control_file_pattern <- paste0(clade, "_NON_SS_extracted_MAFFTA_PAS_context.fa")
  terminal_file_pattern <- paste0(clade, "_SS_extracted_MAFFTA_PAS_context.fa")
  
  # List files in each directory that match the pattern
  control_files <- list.files(control_dir, pattern = control_file_pattern, full.names = TRUE)
  terminal_files <- list.files(terminal_dir, pattern = terminal_file_pattern, full.names = TRUE)
  
  # Read the FASTA files
  control_sequences <- lapply(control_files, function(file) {
    readDNAStringSet(file)
  })
  
  terminal_sequences <- lapply(terminal_files, function(file) {
    readDNAStringSet(file)
  })
  
  # Convert DNAStringSet to character vectors
  control_sequences <- lapply(control_sequences, as.character)
  terminal_sequences <- lapply(terminal_sequences, as.character)
  
  # Flatten the list of character vectors into a single vector for each
  control_sequences <- unlist(control_sequences)
  terminal_sequences <- unlist(terminal_sequences)
  
  return(list(control = control_sequences, terminal = terminal_sequences))
}

# Function to trim sequences
trim_sequences <- function(sequences, trim_start, trim_end) {
  sapply(sequences, function(seq) {
    if (nchar(seq) > (trim_start + trim_end)) {
      substr(seq, trim_start + 1, nchar(seq) - trim_end)
    } else {
      ""  # Return an empty string if the sequence is too short
    }
  })
}

# List of clade names
#clade_names <- c("MLT1A", "MLT1A0", "MLT1A1", "MLT1B", "MLT1C", "MLT1C2")
clade_names <- "MLT1A01BC2"

# Initialize vectors to store all sequences
all_control_sequences <- character()
all_terminal_sequences <- character()

# Loop over each clade name
for (clade_name in clade_names) {
  sequences <- read_fasta_files(clade_name)
  
  # Adjust sequences for specific clades
  #sequences$control <- adjust_sequences(sequences$control, clade_name)
  #sequences$terminal <- adjust_sequences(sequences$terminal, clade_name)
  
  # Trim sequences
  if (clade_name == "MLT1B" || clade_name == "MLT1C") {
    # Trim 15 bases from the start and 15 from the end for MLT1B and MLT1C
    trimmed_control_sequences <- trim_sequences(sequences$control, 15, 14)
    trimmed_terminal_sequences <- trim_sequences(sequences$terminal, 15, 14)
  } else if (clade_name == "MLT1C2") {
    # Trim 15 bases from the start and 13 from the end for MLT1C2
    trimmed_control_sequences <- trim_sequences(sequences$control, 15, 13)
    trimmed_terminal_sequences <- trim_sequences(sequences$terminal, 15, 13)
  } else {
    # Trim 15 bases from the start and end for other clades
    trimmed_control_sequences <- trim_sequences(sequences$control, 15, 15)
    trimmed_terminal_sequences <- trim_sequences(sequences$terminal, 15, 15)
  }
  
  # Extend the vectors with new sequences
  all_control_sequences <- c(all_control_sequences, trimmed_control_sequences)
  all_terminal_sequences <- c(all_terminal_sequences, trimmed_terminal_sequences)
}

# Set the y-axis limit to 2 bits
y_axis_limit <- 2

# Set the y-axis limit to 2 bits
y_axis_limit <- 2

# Create a descriptive title
clade_description <- paste(clade_names, collapse = ", ")
title_control <- "Control Sequences"
title_terminal <- "Terminal Sequences"

# Plot the control sequences
if (length(all_control_sequences) > 0) {
  ggseqlogo(all_control_sequences, method = "bits") +
    ggtitle(title_control) +
    ylim(0, y_axis_limit) +
    scale_x_continuous(breaks = seq(0, max(nchar(all_control_sequences)), by = 5)) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 28, color = "black",hjust = 0.5),
      axis.title = element_text(size = 24, color = "black"),  # Set color to black
      axis.text = element_text(size = 22, color = "black")   # Set color to black
    )
}

# Plot the terminal sequences
if (length(all_terminal_sequences) > 0) {
  ggseqlogo(all_terminal_sequences, method = "bits") +
    ggtitle(title_terminal) +
    ylim(0, y_axis_limit) +
    scale_x_continuous(breaks = seq(0, max(nchar(all_terminal_sequences)), by = 5)) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 28, color = "black",hjust = 0.5),
      axis.title = element_text(size = 24, color = "black"),  # Set color to black
      axis.text = element_text(size = 22, color = "black")   # Set color to black
    )
}

### STATISTICAL TESTING ###
# Function to calculate Shannon entropy for a set of sequences
calculate_entropy <- function(sequences) {
  # Convert sequences to a matrix
  seq_matrix <- do.call(rbind, strsplit(sequences, ""))
  
  # Calculate frequency of each nucleotide at each position
  freq_matrix <- apply(seq_matrix, 2, function(column) {
    table(factor(column, levels = c("A", "C", "G", "T"))) / length(column)
  })
  
  # Calculate entropy for each position
  entropy <- apply(freq_matrix, 2, function(freqs) {
    -sum(freqs * log2(freqs + 1e-10))  # Add a small value to avoid log(0)
  })
  
  return(entropy)
}
# Function to trim sequences to a specific length
trim_for_test <- function(sequences, trim_start, trim_end) {
  sapply(sequences, function(seq) {
    if (nchar(seq) > (trim_start + trim_end)) {
      substr(seq, trim_start + 1, nchar(seq) - trim_end)
    } else {
      ""  # Return an empty string if the sequence is too short
    }
  })
}


# Plot entropy profiles
plot(control_entropy, type = "l", col = "blue", ylim = c(0, 2), ylab = "Entropy", xlab = "Position", main = "Shannon Entropy Profiles")
lines(terminal_entropy, col = "red")
legend("topright", legend = c("Control", "Terminal"), col = c("blue", "red"), lty = 1)

# Trim sequences for the Wilcoxon test
trimmed_control_for_test <- trim_for_test(all_control_sequences, 15, 15)
trimmed_terminal_for_test <- trim_for_test(all_terminal_sequences, 15, 15)

# Filter out empty sequences that may result from trimming
trimmed_control_for_test <- trimmed_control_for_test[nchar(trimmed_control_for_test) == 6]
trimmed_terminal_for_test <- trimmed_terminal_for_test[nchar(trimmed_terminal_for_test) == 6]

# Calculate entropy for trimmed control and terminal sequences
control_entropy <- calculate_entropy(trimmed_control_for_test)
terminal_entropy <- calculate_entropy(trimmed_terminal_for_test)


# Function to trim sequences to a specific length
trim_for_test <- function(sequences, trim_start, trim_end) {
  sapply(sequences, function(seq) {
    if (nchar(seq) > (trim_start + trim_end)) {
      substr(seq, trim_start + 1, nchar(seq) - trim_end)
    } else {
      ""  # Return an empty string if the sequence is too short
    }
  })
}

# Trim sequences for the Wilcoxon test
trimmed_control_for_test <- trim_for_test(all_control_sequences, 15, 15)
trimmed_terminal_for_test <- trim_for_test(all_terminal_sequences, 15, 15)

# Filter out empty sequences that may result from trimming
trimmed_control_for_test <- trimmed_control_for_test[nchar(trimmed_control_for_test) == 6]
trimmed_terminal_for_test <- trimmed_terminal_for_test[nchar(trimmed_terminal_for_test) == 6]

# Calculate entropy for trimmed control and terminal sequences
control_entropy <- calculate_entropy(trimmed_control_for_test)
terminal_entropy <- calculate_entropy(trimmed_terminal_for_test)

# Perform a Wilcoxon signed-rank test to compare entropy
test_result <- wilcox.test(control_entropy, terminal_entropy, paired = TRUE)

# Annotate the plot with the Wilcoxon test results in the bottom left
text(x = 1, y = 0.1, 
     labels = paste("Wilcoxon p-value:", format(test_result$p.value, digits = 3)), 
     col = "black", cex = 0.8, adj = c(0, 0))