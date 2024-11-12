#!/bin/bash
DATABASES_DIR="inputs/databases/"
REF_SPECIES_PATH="inputs/reference/"
REF_SPECIES_NO_EXT_PATH="inputs/reference/"
NONREF_SPECIES_DIR="/home/amohs002/projects/research/ALLEGRO_Fungi_Downloader/data/fourdbs_concat/proteomes/"
OUTPUT_DIR="outputs/"

# Find the file in the directory that does not have the .dmnd extension
for file in "$REF_SPECIES_PATH"/*; do
    if [[ -f "$file" && "$file" != *.dmnd ]]; then
        # Extract the filename
        filename_with_ext=$(basename "$file")
        REF_SPECIES_PATH=$REF_SPECIES_PATH$filename_with_ext
        REF_SPECIES_NO_EXT_PATH=$REF_SPECIES_NO_EXT_PATH${filename_with_ext%.*}
        break
    fi
done

echo "Query Non-Ref vs Ref"
for FILE in $(ls $NONREF_SPECIES_DIR)
do
    SPECIES_NAME="$(echo "$FILE" | cut -d'.' -f1)"
    OUT=${OUTPUT_DIR}${SPECIES_NAME}_vs_ref.tsv

    diamond blastp -q $NONREF_SPECIES_DIR$FILE -d $REF_SPECIES_NO_EXT_PATH -o $OUT --ultra-sensitive --top 0 --quiet
done

echo "Query Ref vs Non-Ref"
for DB in $(ls $DATABASES_DIR)
do
    SPECIES_NAME="$(echo "$DB" | cut -d'.' -f1)"
    DB_NAME="${DATABASES_DIR}${SPECIES_NAME}"
    OUT="${OUTPUT_DIR}ref_vs_${SPECIES_NAME}.tsv"

    diamond blastp -q $REF_SPECIES_PATH -d $DB_NAME -o $OUT --ultra-sensitive --top 0 --quiet
done
echo "Done."
