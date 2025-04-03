from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import os
import re

def process_fasta(input_file, output_file):
    name_count = defaultdict(int)
    with open(output_file, "w") as out_f:
        for record in SeqIO.parse(input_file, "fasta"):
            header = record.id
            chrom_band, herv_clade, sense = None, None, None
            
            # Use regular expressions to parse the header
            match = re.search(r"-(\d+)_([^_]+(?:_[^_]+)*)_(\+|-)", header)
            if match:
                chrom_band = header.split("_")[0]  # Extract the chrom_band
                herv_clade = match.group(2)
                sense = match.group(3)
            
            # Construct the new header with only cytoband and HERV clade
            new_header = f"{chrom_band}_{herv_clade}"
            if sense == "-":
                new_header += "_RC"
                record.seq = record.seq.reverse_complement()
            
            # Ensure unique headers
            name_count[new_header] += 1
            if name_count[new_header] > 1:
                new_header += f"_{name_count[new_header]}"
                
            record.id = new_header
            record.description = ""
            SeqIO.write(record, out_f, "fasta")

# Define the arrays for SS_or_NON_SS and Clades
SS_or_NON_SS_array = ["SS", "NON_SS"]

# Read the clades from the CladeList.tsv file
clade_list_file = "/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs/Last_Exon_Clades/CladeList.tsv"
with open(clade_list_file, "r") as f:
    Clades = [line.strip() for line in f]

# Base directory for input and output files
base_dir = "/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs"

for SS_or_NON_SS in SS_or_NON_SS_array:
    for Clade in Clades:
        if SS_or_NON_SS == "SS":
            input_dir = f"{base_dir}/terminal_HERV_fastas"
            output_dir = f"{base_dir}/terminal_HERV_fastas"
        else:
            input_dir = f"{base_dir}/Control_NON_terminal_HERV_fastas"
            output_dir = f"{base_dir}/Control_NON_terminal_HERV_fastas"
        
        input_file = f"{input_dir}/3prmTerminalHERVs/3prmTerminalHERVseqs_{Clade}.fa"
        output_file = f"{output_dir}/First_Step_Processed/processed_{Clade}_{SS_or_NON_SS}_3prmTerminalHERVseqs.fa"
        
        # Check if the input file exists
        if os.path.exists(input_file):
            print(f"Processing {input_file} to {output_file}")
            process_fasta(input_file, output_file)
        else:
            print(f"Input file not found: {input_file}")