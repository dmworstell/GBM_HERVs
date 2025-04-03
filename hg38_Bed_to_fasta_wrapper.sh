#!/bin/bash

Genome=hg38

#conda activate bedtools_env

SRCE_dir=/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project

Genome_fasta_dir=${SRCE_dir}/Reference/Genomes

Genome_fasta=$Genome_fasta_dir/hg38.p14.fa

IN_dir=${SRCE_dir}/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs

SS_or_NON_SS_array=("SS" "NON_SS")

# Read the clades from the CladeList.tsv file
CladeList_file=${IN_dir}/Last_Exon_Clades/CladeList.tsv
Clades=($(cut -f1 $CladeList_file))

# Create a temporary file to store valid clades
Valid_Clades_file=$(mktemp)

for Clade in "${Clades[@]}"; do
    SS_bed=${IN_dir}/Last_Exon_Clades/DE_SS_LE_PAS_${Clade}.bed
    NON_SS_bed=${IN_dir}/Last_Exon_Clades/DE_NON_SS_LE_PAS_${Clade}.bed

    # Check if both SS and NON_SS bed files exist
    if [[ -f "$SS_bed" && -f "$NON_SS_bed" ]]; then
        echo "$Clade" >> "$Valid_Clades_file"
    else
        echo "Skipping clade ${Clade} as one or both bed files are missing."
    fi
done

# Replace the original CladeList.tsv with the valid clades
mv "$Valid_Clades_file" "$CladeList_file"

# Read the updated clades from the CladeList.tsv file
Clades=($(cut -f1 $CladeList_file))

for SS_or_NON_SS in "${SS_or_NON_SS_array[@]}"; do
    for Clade in "${Clades[@]}"; do
        # Set OUT_dir based on the value of SS_or_NON_SS
        if [ "$SS_or_NON_SS" = "SS" ]; then
            OUT_dir=${IN_dir}/terminal_HERV_fastas
        else
            OUT_dir=${IN_dir}/Control_NON_terminal_HERV_fastas
        fi

        IN_bed=${IN_dir}/Last_Exon_Clades/DE_${SS_or_NON_SS}_LE_PAS_${Clade}.bed

        echo "Processing ${SS_or_NON_SS} with Clade ${Clade}"
        echo BEGIN

        bedtools getfasta -fi ${Genome_fasta} -bed ${IN_bed} -fo ${OUT_dir}/3prmterminalHERVs/3prmTerminalHERVseqs_${Clade}.fa -name

        echo DONE
    done
done

#conda deactivate
