#!/bin/bash

# Define the directory containing the .bed and .tsv files
TARGET_DIR="/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs/Last_Exon_Clades"
# Define the path to the cladelist file
CLADELIST_FILE="/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs/Last_Exon_Clades/CladeList.tsv"

# Loop through all .bed and .tsv files in the target directory
for FILE in "$TARGET_DIR"/*.{bed,tsv}; do
    # Check if the file name contains a clade starting with HERV or ERV and ending with -int
    if [[ "$FILE" =~ (HERV|ERV).*[-_]int\.(bed|tsv)$ ]]; then
        # Extract the directory, filename, and extension
        DIR=$(dirname "$FILE")
        FILENAME=$(basename "$FILE")
        EXTENSION="${FILENAME##*.}"
        
        # Remove the -int suffix from the filename
        NEW_FILENAME=$(echo "$FILENAME" | sed 's/-int//')
        NEW_FILE="$DIR/$NEW_FILENAME"
        
        # Rename the file, overwriting if necessary
        mv -f "$FILE" "$NEW_FILE"
        echo "Renamed $FILE to $NEW_FILE"
        
        # Update the contents of the file
        sed -E -i '' 's/(HUERS[a-zA-Z0-9_-]*|HERV[a-zA-Z0-9_-]*|ERV[a-zA-Z0-9_-]*)-int/\1/g' "$NEW_FILE"
        echo "Updated contents of $NEW_FILE"
    fi
done

# Update the cladelist file
if [ -f "$CLADELIST_FILE" ]; then
    sed -E -i '' 's/(HUERS[a-zA-Z0-9_-]*|HERV[a-zA-Z0-9_-]*|ERV[a-zA-Z0-9_-]*)-int/\1/g' "$CLADELIST_FILE"
    echo "Updated cladelist file: $CLADELIST_FILE"
else
    echo "Cladelist file not found: $CLADELIST_FILE"
fi

echo "All files processed successfully."
