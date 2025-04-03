#!/bin/bash

# Define the base directory and clade list file
BASE_DIR="/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs"
CLADE_LIST_FILE="${BASE_DIR}/Last_Exon_Clades/CladeList.tsv"

# Define the directories for input and output files
COMBINED_ALIGNMENTS_DIR="${BASE_DIR}/combined_alignments"
TERMINAL_HERV_DIR="${BASE_DIR}/terminal_HERV_fastas/First_Step_Processed"
CONTROL_NON_TERMINAL_HERV_DIR="${BASE_DIR}/Control_NON_terminal_HERV_fastas/First_Step_Processed"
SS_OUTPUT_DIR="${BASE_DIR}/terminal_HERV_fastas/Aligned_HERVs"
NON_SS_OUTPUT_DIR="${BASE_DIR}/Control_NON_terminal_HERV_fastas/Aligned_HERVs"

# Create the output directories if they don't exist
mkdir -p "$SS_OUTPUT_DIR"
mkdir -p "$NON_SS_OUTPUT_DIR"

# Read the clades from the CladeList.tsv file
Clades=($(cut -f1 $CLADE_LIST_FILE))

# Function to extract headers from a fasta file
extract_headers() {
    grep "^>" "$1" | sed 's/^>//'
}

# Loop through each clade and process the combined alignment file
for Clade in "${Clades[@]}"; do
    # Define the input files for SS and NON_SS
    SS_FILE="${TERMINAL_HERV_DIR}/processed_${Clade}_SS_3prmTerminalHERVseqs.fa"
    NON_SS_FILE="${CONTROL_NON_TERMINAL_HERV_DIR}/processed_${Clade}_NON_SS_3prmTerminalHERVseqs.fa"
    
    # Define the combined alignment file
    COMBINED_FILE="${COMBINED_ALIGNMENTS_DIR}/combined_alignment_${Clade}.fa"
    
    # Define the output files
    SS_OUTPUT_FILE="${SS_OUTPUT_DIR}/aligned_${Clade}_SS_3prmTerminalHERVseqs.fa"
    NON_SS_OUTPUT_FILE="${NON_SS_OUTPUT_DIR}/aligned_${Clade}_NON_SS_3prmTerminalHERVseqs.fa"
    
    # Extract headers from the original SS and NON_SS files
    SS_HEADERS=$(extract_headers "$SS_FILE")
    NON_SS_HEADERS=$(extract_headers "$NON_SS_FILE")
    
    # Convert headers to arrays
    IFS=$'\n' read -d '' -r -a ss_headers_array <<< "$SS_HEADERS"
    IFS=$'\n' read -d '' -r -a non_ss_headers_array <<< "$NON_SS_HEADERS"
    
    # Initialize output files
    > "$SS_OUTPUT_FILE"
    > "$NON_SS_OUTPUT_FILE"
    
    # Read the combined alignment file and process sequences
    awk -v ss_output="$SS_OUTPUT_FILE" -v non_ss_output="$NON_SS_OUTPUT_FILE" '
    BEGIN {
        for (i in ARGV) {
            if (ARGV[i] ~ /^ss_header_/) {
                header = substr(ARGV[i], 11)  # Remove the prefix "ss_header_"
                ss_map[header] = 1
                delete ARGV[i]
            } else if (ARGV[i] ~ /^non_ss_header_/) {
                header = substr(ARGV[i], 15)  # Remove the prefix "non_ss_header_"
                non_ss_map[header] = 1
                delete ARGV[i]
            }
        }
    }
    /^>/ {
        header = substr($0, 2)
        split(header, parts, "|")
        stripped_header = parts[2]
        if (stripped_header in ss_map) {
            output = ss_output
        } else if (stripped_header in non_ss_map) {
            output = non_ss_output
        } else {
            output = ""
        }
    }
    {
        if (output != "") {
            print $0 > output
        }
    }
    ' "${ss_headers_array[@]/#/ss_header_}" "${non_ss_headers_array[@]/#/non_ss_header_}" "$COMBINED_FILE"
    
    echo "Processed clade ${Clade}. SS output saved to ${SS_OUTPUT_FILE}, NON_SS output saved to ${NON_SS_OUTPUT_FILE}."
done

echo "All clades processed successfully."
