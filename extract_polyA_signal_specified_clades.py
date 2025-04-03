#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from Bio import SeqIO

# Define the base directory
BASE_DIR = "/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs"

# Define the directories for input and output files
TERMINAL_HERV_DIR = os.path.join(BASE_DIR, "terminal_HERV_fastas/Aligned_HERVs/specified_clades")
CONTROL_NON_TERMINAL_HERV_DIR = os.path.join(BASE_DIR, "Control_NON_terminal_HERV_fastas/Aligned_HERVs/specified_clades")
CONSENSUS_DIR = os.path.join(BASE_DIR, "HERV_consensus_sequences/extracted_consensus_sequences")
SS_OUTPUT_DIR = os.path.join(BASE_DIR, "terminal_HERV_fastas/polyA_contexts/specified_clades")
NON_SS_OUTPUT_DIR = os.path.join(BASE_DIR, "Control_NON_terminal_HERV_fastas/polyA_Contexts/specified_clades")

# Create the output directories if they don't exist
os.makedirs(SS_OUTPUT_DIR, exist_ok=True)
os.makedirs(NON_SS_OUTPUT_DIR, exist_ok=True)

# Define the combined clades alignment files
SS_FILE = os.path.join(TERMINAL_HERV_DIR, "aligned_MLT1A01BC2_SS.fa")
NON_SS_FILE = os.path.join(CONTROL_NON_TERMINAL_HERV_DIR, "aligned_MLT1A01BC2_NON_SS.fa")

# Define the output files
SS_OUTPUT_FILE = os.path.join(SS_OUTPUT_DIR, "MLT1A01BC2_SS_extracted_MAFFTA_PAS_context.fa")
NON_SS_OUTPUT_FILE = os.path.join(NON_SS_OUTPUT_DIR, "MLT1A01BC2_NON_SS_extracted_MAFFTA_PAS_context.fa")

# Define the consensus sequence file for MLT1A
CONSENSUS_FILE = os.path.join(CONSENSUS_DIR, "MLT1A.fasta")

# Function to extract polyA signal context
def extract_polyA_signal_context(input_file, output_file, start, end):
    with open(output_file, "w") as out_f:
        for record in SeqIO.parse(input_file, "fasta"):
            # Extract the sequence from positions start to end (1-based indexing)
            extracted_seq = record.seq[start-1:end]
            # Pad the sequence if it's too short
            if len(extracted_seq) < 66:
                if start <= 30:
                    extracted_seq = '-' * (30 - start + 1) + extracted_seq
                if end >= len(record.seq) - 30:
                    extracted_seq = extracted_seq + '-' * (66 - len(extracted_seq))
            # Ensure the extracted sequence is exactly 66 bases long
            extracted_seq = extracted_seq[:66]
            # Create a new record with the extracted sequence
            new_record = record[:0]  # Create an empty record with the same ID
            new_record.seq = extracted_seq
            new_record.id = record.id
            new_record.description = ""
            # Write the new record to the output file
            SeqIO.write(new_record, out_f, "fasta")

# Function to find the polyA signal in the consensus sequence
def find_polyA_signal(consensus_seq):
    canonical_signal = "AATAAA"
    half_length = len(consensus_seq) // 2
    positions = [i for i in range(half_length, len(consensus_seq) - 5) if consensus_seq[i:i+6].upper() == canonical_signal]

    if positions:
        return positions[-1], canonical_signal  # Return the most downstream position in the last half

    return None, None

# Read the MLT1A consensus sequence
consensus_seq = str(next(SeqIO.parse(CONSENSUS_FILE, "fasta")).seq)

# Find the polyA signal in the MLT1A consensus sequence
pos, signal = find_polyA_signal(consensus_seq)

if pos is not None:
    start = max(1, pos - 30 + 1)  # Ensure start is at least 1
    end = min(len(consensus_seq), pos + 6 + 30)  # Corrected end position calculation
    # Extract polyA signal context for SS and NON_SS files
    extract_polyA_signal_context(SS_FILE, SS_OUTPUT_FILE, start, end)
    extract_polyA_signal_context(NON_SS_FILE, NON_SS_OUTPUT_FILE, start, end)
    print(f"Processed combined clades. SS output saved to {SS_OUTPUT_FILE}, NON_SS output saved to {NON_SS_OUTPUT_FILE}.")
else:
    print("No polyA signal found in the last half for MLT1A consensus sequence. Skipping.")

print("All specified clades processed successfully.")
