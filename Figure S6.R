library(ape)
library(phangorn)
library(msa)

# Define the directory containing the .fasta files
directory <- "/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs/HERV_consensus_sequences/extracted_consensus_sequences"

# List all .fasta files
fasta_files <- list.files(directory, pattern = "\\.fasta$", full.names = TRUE)

# Read all sequences into a list and assign names based on file names
sequences <- lapply(fasta_files, function(file) {
  # Read the sequence as a character vector
  seq <- readLines(file)
  
  # Remove any lines that are headers (start with '>')
  seq <- seq[!grepl("^>", seq)]
  
  # Concatenate the sequence lines and remove spaces
  seq <- gsub(" ", "", paste(seq, collapse = ""))
  
  return(seq)
})

# Extract clade names from filenames
clade_names <- gsub("\\.fasta$", "", basename(fasta_files))

# Combine all sequences into a single object
all_sequences <- do.call(c, sequences)
names(all_sequences) <- clade_names

# Subset sequences to only those starting with "MLT1"
subset_indices <- grep("^(MLT1)", clade_names)
subset_sequences <- all_sequences[subset_indices]

# Debug: Print the names of the subset sequences
cat("Subset sequence names:\n")
print(names(subset_sequences))

# Check if there are enough sequences to align
if (length(subset_sequences) < 2) {
  stop("Not enough sequences to perform alignment.")
}

# Perform alignment using msa (which uses MUSCLE internally)
alignment <- msa(subset_sequences, method = "Muscle", type = "dna")

# Print the alignment to check
print(alignment)

# Convert alignment to a format suitable for distance calculation
alignment_dnabin <- as.DNAbin(alignment)

# Compute a distance matrix
dist_matrix <- dist.dna(alignment_dnabin, model = "raw")

# Convert the distance object to a full matrix
dist_matrix_full <- as.matrix(dist_matrix)

# Build a phylogenetic tree using Neighbor-Joining
tree <- nj(dist_matrix)

# Create a vector of colors for the tip labels
tip_colors <- rep("black", length(tree$tip.label))
names(tip_colors) <- tree$tip.label

# Assign colors to specific clades
tip_colors[c("MLT1A", "MLT1A0", "MLT1A1")] <- "red"
tip_colors["MLT1B"] <- "blue"
tip_colors[c("MLT1C", "MLT1C2")] <- "green"

# Plot the tree with a specific title and colored tip labels
plot(tree, main = "Phylogenetic Tree of MLT1 Consensus Sequences", tip.color = tip_colors)

# Add a scale bar to the plot
add.scale.bar(x = max(tree$edge.length) - 0.05, y = 1.5, length = 0.05, lwd = 2, col = "black")

# Add a label for the scale bar
text(x = max(tree$edge.length) - 0.06, y = 2.5, labels = "Sequence Distance", pos = 4)

# Save the tree to a file
write.tree(tree, file = "phylogenetic_tree.nwk")

### FIND CLOSE TO CLADE ###

# Identify clades close to "MLT1A"
target_clade <- "MLT1A"
threshold_distance <- 0.2  # Define a threshold distance

# Find the index of the target clade
target_index <- which(names(subset_sequences) == target_clade)

# Extract distances to the target clade from the full matrix
distances_to_target <- dist_matrix_full[target_index, ]

# Identify clades within the threshold distance
close_clades <- names(subset_sequences)[distances_to_target <= threshold_distance]

# Print the close clades
cat("Clades close to", target_clade, "within distance", threshold_distance, ":\n")
print(close_clades)