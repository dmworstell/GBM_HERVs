#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
from Bio import SeqIO

# Define the base directory and clade list file
BASE_DIR = "/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs"
SS_INPUT_DIR = os.path.join(BASE_DIR, "terminal_HERV_fastas/polyA_contexts")
NON_SS_INPUT_DIR = os.path.join(BASE_DIR, "Control_NON_terminal_HERV_fastas/polyA_Contexts")
SS_OUTPUT_DIR = os.path.join(BASE_DIR, "terminal_HERV_fastas/combined_polyA_contexts")
NON_SS_OUTPUT_DIR = os.path.join(BASE_DIR, "Control_NON_terminal_HERV_fastas/combined_polyA_contexts")

# Create the output directories if they don't exist
os.makedirs(SS_OUTPUT_DIR, exist_ok=True)
os.makedirs(NON_SS_OUTPUT_DIR, exist_ok=True)

# Function to combine polyA signal context files
def combine_polyA_signal_contexts(prefix, ss_output_file, non_ss_output_file):
    with open(ss_output_file, "w") as ss_out_f, open(non_ss_output_file, "w") as non_ss_out_f:
        for filename in os.listdir(SS_INPUT_DIR):
            if (not prefix or filename.startswith(prefix)) and filename.endswith("_SS_extracted_MAFFTA_PAS_context.fa"):
                with open(os.path.join(SS_INPUT_DIR, filename), "r") as ss_in_f:
                    lines = ss_in_f.readlines()
                    if lines:  # Check if the file is not empty
                        ss_out_f.writelines(lines)
        
        for filename in os.listdir(NON_SS_INPUT_DIR):
            if (not prefix or filename.startswith(prefix)) and filename.endswith("_NON_SS_extracted_MAFFTA_PAS_context.fa"):
                with open(os.path.join(NON_SS_INPUT_DIR, filename), "r") as non_ss_in_f:
                    lines = non_ss_in_f.readlines()
                    if lines:  # Check if the file is not empty
                        non_ss_out_f.writelines(lines)

# Main function to parse arguments and call the combine function
def main():
    parser = argparse.ArgumentParser(description="Combine polyA signal context files based on clade prefix.")
    parser.add_argument("prefix", nargs='?', default="", help="Prefix to filter clades (e.g., MLT1, THE1, LTR). If not supplied, all clades will be combined.")
    args = parser.parse_args()

    prefix = args.prefix
    ss_output_file = os.path.join(SS_OUTPUT_DIR, f"combined_polyA_signal_{prefix}_SS.fa" if prefix else "combined_polyA_signal_SS.fa")
    non_ss_output_file = os.path.join(NON_SS_OUTPUT_DIR, f"combined_polyA_signal_{prefix}_NON_SS.fa" if prefix else "combined_polyA_signal_NON_SS.fa")

    combine_polyA_signal_contexts(prefix, ss_output_file, non_ss_output_file)
    print(f"Combined polyA signal context files for prefix '{prefix}' saved to {ss_output_file} and {non_ss_output_file}.")

if __name__ == "__main__":
    main()