#!/bin/bash

# Define the base directory and clade list file
BASE_DIR="/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs"
CLADE_LIST_FILE="${BASE_DIR}/Last_Exon_Clades/CladeList.tsv"

# Define the directories for input and output files
TERMINAL_HERV_DIR="${BASE_DIR}/terminal_HERV_fastas/First_Step_Processed"
CONTROL_NON_TERMINAL_HERV_DIR="${BASE_DIR}/Control_NON_terminal_HERV_fastas/First_Step_Processed"
OUTPUT_DIR="${BASE_DIR}/combined_alignments"
CONSENSUS_LIB_FILE="${BASE_DIR}/HERV_consensus_sequences/Dfam-RepeatMasker.lib"
EXTRACTED_CONSENSUS_DIR="${BASE_DIR}/HERV_consensus_sequences/extracted_consensus_sequences"

# Create the output directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$EXTRACTED_CONSENSUS_DIR"

# Read the clades from the CladeList.tsv file
Clades=($(cut -f1 $CLADE_LIST_FILE))

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

# Loop through each clade and run the mafft command
for Clade in "${Clades[@]}"; do
    # Define the input files for SS and NON_SS
    SS_FILE="${TERMINAL_HERV_DIR}/processed_${Clade}_SS_3prmTerminalHERVseqs.fa"
    NON_SS_FILE="${CONTROL_NON_TERMINAL_HERV_DIR}/processed_${Clade}_NON_SS_3prmTerminalHERVseqs.fa"
    
    # Define the temporary consensus sequence file for the clade
    CONSENSUS_FILE="${EXTRACTED_CONSENSUS_DIR}/${Clade}.fasta"
    
    # Extract the consensus sequence from the library file
    extract_consensus_sequence "$Clade" "$CONSENSUS_FILE"
    
    # Define the output file
    OUTPUT_FILE="${OUTPUT_DIR}/combined_alignment_${Clade}.fa"
    
    # Create a temporary file for the combined input
    TEMP_INPUT_FILE=$(mktemp)
    
    # Combine the input files into the temporary file
    cat "$SS_FILE" "$NON_SS_FILE" > "$TEMP_INPUT_FILE"
    
    # Print the path of the temporary file for inspection
    echo "Temporary input file for clade ${Clade}: $TEMP_INPUT_FILE"
    # Combine the input files and run the mafft command
    echo "Running mafft for clade ${Clade}..."
    mafft --auto --reorder --keeplength --compactmapout --anysymbol --addfragments "$TEMP_INPUT_FILE" "$CONSENSUS_FILE" > "$OUTPUT_FILE"
    
    if [ $? -ne 0 ]; then
        echo "Error: mafft command failed for clade ${Clade}."
        exit 1
    fi
    
    echo "Completed mafft for clade ${Clade}. Output saved to ${OUTPUT_FILE}."
    # Remove the temporary file
    rm "$TEMP_INPUT_FILE"
done

echo "All clades processed successfully."
