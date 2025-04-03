#!/bin/bash

# Define the directories to check
DIR1="/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs/terminal_HERV_fastas/polyA_contexts"
DIR2="/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs/Control_NON_terminal_HERV_fastas/polyA_Contexts"

# Function to check fasta files in a directory
check_fasta_files() {
    local DIR=$1
    echo "Checking directory: $DIR"
    for FILE in "$DIR"/*.fa; do
        if grep -q -E "NO POLYA SIGNAL FOUND|MULTIPLE" "$FILE"; then
            # Extract the clade name from the file path
            FILENAME=$(basename "$FILE")
            CLADE=$(echo "$FILENAME" | cut -d'_' -f1)
            echo "Found in clade: $CLADE"
        fi
    done
}

# Check the directories
check_fasta_files "$DIR1"
check_fasta_files "$DIR2"

echo "Check completed."
