#!/bin/bash

# Define the directory containing the scripts
SCRIPT_DIR="/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/Scripts/polyA_site_scripts"

# Run the hg38_Bed_to_fasta_wrapper.sh script
echo "Running hg38_Bed_to_fasta_wrapper.sh..."
bash "$SCRIPT_DIR/hg38_Bed_to_fasta_wrapper.sh"
if [ $? -ne 0 ]; then
    echo "Error: hg38_Bed_to_fasta_wrapper.sh failed."
    exit 1
fi

# Run the Manipulate_Fasta.py script
echo "Running Manipulate_Fasta.py..."
python3 "$SCRIPT_DIR/Manipulate_Fasta.py"
if [ $? -ne 0 ]; then
    echo "Error: Manipulate_Fasta.py failed."
    exit 1
fi

# Run the run_mafft.sh script
echo "Running run_mafft.sh..."
bash "$SCRIPT_DIR/run_mafft.sh"
if [ $? -ne 0 ]; then
    echo "Error: run_mafft.sh failed."
    exit 1
fi

# Run the split_alignments.sh script
echo "Running split_alignments.sh..."
bash "$SCRIPT_DIR/split_alignments.sh"
if [ $? -ne 0 ]; then
    echo "Error: split_alignments.sh failed."
    exit 1
fi

# Run the extract_polyA_signal.py script
echo "Running extract_polyA_signal.py..."
python3 "$SCRIPT_DIR/extract_polyA_signal.py"
if [ $? -ne 0 ]; then
    echo "Error: extract_polyA_signal.py failed."
    exit 1
fi

echo "All scripts ran successfully."

# Run the combine_polyA_signal_contexts.py script
echo "Running combine_polyA_signal_contexts.py..."
python3 "$SCRIPT_DIR/combine_polyA_signal_contexts.py"
if [ $? -ne 0 ]; then
    echo "Error: combine_polyA_signal_contexts.py"
    exit 1
fi

echo "All scripts ran successfully."
