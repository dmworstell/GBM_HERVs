#!/bin/bash

# Define the directory containing the .bed and .tsv files
TARGET_DIR="/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs/Last_Exon_Clades"
# Define the path to the CladeList.tsv file
CLADELIST_FILE="$TARGET_DIR/CladeList.tsv"

# Loop through all .bed and .tsv files in the target directory
for FILE in "$TARGET_DIR"/*.{bed,tsv}; do
    # Check if the file name ends with -int
    if [[ "$FILE" =~ (-int|_I|-I)\.(bed|tsv)$ ]]; then
        # Remove the file
        rm "$FILE"
        echo "Removed $FILE"
    fi
done

# Update the CladeList.tsv file
if [ -f "$CLADELIST_FILE" ]; then
    sed -E -i '' '/-int/d' "$CLADELIST_FILE"
    sed -E -i '' '/_I/d' "$CLADELIST_FILE"
    sed -E -i '' '/-I/d' "$CLADELIST_FILE"
    echo "Updated CladeList.tsv file: $CLADELIST_FILE"
else
    echo "CladeList.tsv file not found: $CLADELIST_FILE"
fi

echo "All files processed successfully."
