#!/bin/bash

# Define the base directory and clade list file
BASE_DIR="/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs"

# Define the directories for input and output files
TERMINAL_HERV_DIR="${BASE_DIR}/terminal_HERV_fastas/First_Step_Processed"
CONTROL_NON_TERMINAL_HERV_DIR="${BASE_DIR}/Control_NON_terminal_HERV_fastas/First_Step_Processed"
OUTPUT_DIR="${BASE_DIR}/combined_alignments/specified_clades"
CONSENSUS_LIB_FILE="${BASE_DIR}/HERV_consensus_sequences/Dfam-RepeatMasker.lib"
EXTRACTED_CONSENSUS_DIR="${BASE_DIR}/HERV_consensus_sequences/extracted_consensus_sequences"

# Create the output and consensus directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$EXTRACTED_CONSENSUS_DIR"

# Define the subset of clades to process
SelectedClades=("MLT1A" "MLT1A0" "MLT1B" "MLT1C" "MLT1C2")

# Define the output file
OUTPUT_FILE="${OUTPUT_DIR}/combined_alignment_MLT1A01BC2.fa"

# Create a temporary file for the combined input
TEMP_INPUT_FILE=$(mktemp)

# Function to extract consensus sequence from the library file
extract_consensus_sequence() {
    local clade=$1
    local consensus_file=$2
    local found=0

    while IFS= read -r line; do
        if [[ $line == ">$clade#LTR"* ]]; then
            found=1
            echo "$line" > "$consensus_file"
        elif [[ $found -eq 1 ]]; then
            if [[ $line == ">"* ]]; then
                break
            fi
            echo "$line" >> "$consensus_file"
        fi
    done < "$CONSENSUS_LIB_FILE"

    if [[ $found -eq 0 ]]; then
        echo "Error: Consensus sequence for clade ${clade} not found in ${CONSENSUS_LIB_FILE}."
        exit 1
    fi
}

# Extract the MLT1A consensus sequence
MLT1A_CONSENSUS_FILE="${EXTRACTED_CONSENSUS_DIR}/MLT1A.fasta"
extract_consensus_sequence "MLT1A" "$MLT1A_CONSENSUS_FILE"

# Loop through each selected clade and combine the input files
for Clade in "${SelectedClades[@]}"; do
    # Define the input files for SS and NON_SS
    SS_FILE="${TERMINAL_HERV_DIR}/processed_${Clade}_SS_3prmTerminalHERVseqs.fa"
    NON_SS_FILE="${CONTROL_NON_TERMINAL_HERV_DIR}/processed_${Clade}_NON_SS_3prmTerminalHERVseqs.fa"
    
    # Check if the files exist before concatenating
    if [[ -f "$SS_FILE" ]]; then
        cat "$SS_FILE" >> "$TEMP_INPUT_FILE"
    else
        echo "Warning: SS file for clade ${Clade} not found."
    fi
    
    if [[ -f "$NON_SS_FILE" ]]; then
        cat "$NON_SS_FILE" >> "$TEMP_INPUT_FILE"
    else
        echo "Warning: NON_SS file for clade ${Clade} not found."
    fi
done

# Run the mafft command with --addfragments
echo "Running mafft with --addfragments on combined input for specified clades..."
mafft --auto --reorder --keeplength --compactmapout --anysymbol --addfragments "$TEMP_INPUT_FILE" "$MLT1A_CONSENSUS_FILE" > "$OUTPUT_FILE"

if [ $? -ne 0 ]; then
    echo "Error: mafft command failed."
    exit 1
fi

echo "Completed mafft. Output saved to ${OUTPUT_FILE}."

# Remove the temporary file
rm "$TEMP_INPUT_FILE"

echo "All specified clades processed successfully."
