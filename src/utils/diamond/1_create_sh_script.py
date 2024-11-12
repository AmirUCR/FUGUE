import os
import re
import yaml
import os
import stat

if not os.path.exists('inputs'):
    os.makedirs('inputs/reference')
    os.makedirs('inputs/databases')

if not os.path.exists('inputs/reference'):
    os.mkdir('inputs/reference')

if not os.path.exists('inputs/databases'):
    os.mkdir('inputs/databases')

if not os.path.exists('outputs'):
    os.mkdir('outputs')

if not os.path.exists('outputs/reciprocal'):
    os.mkdir('outputs/reciprocal')

# ----

yaml_file_path = 'config.yaml'

with open(yaml_file_path, 'r') as file:
    config = yaml.safe_load(file)

proteome_path = config['proteome_path']
match = re.search(r'([^/]+)\.[^/]+$', proteome_path)
if match:
    filename_without_ext = match.group(1)
else:
    filename_without_ext = 'reference_species_db'

sh_content = f"""#!/bin/bash
# Generated script by 1_create_sh_script.py

PROT_PATH="{proteome_path}"
FILENAME={filename_without_ext}

diamond makedb --in $PROT_PATH -d inputs/reference/$FILENAME
"""

sh_file_path = '3_makedb_for_ref_species.sh'

with open(sh_file_path, 'w') as file:
    file.write(sh_content)

# Give exec permissions
st = os.stat(sh_file_path)
os.chmod(sh_file_path, st.st_mode | stat.S_IEXEC)



sh_content = f"""#!/bin/bash
# Generated script by 1_create_sh_script.py
DIR="{config['query_proteins_dir']}"

for FILE in $(ls $DIR)
do
    NAME=$(echo "$FILE" | cut -d'.' -f1)
    diamond makedb --in $DIR$FILE -d inputs/databases/$NAME
done
"""

sh_file_path = '4_makedb_for_all_species.sh'

with open(sh_file_path, 'w') as file:
    file.write(sh_content)

# Give exec permissions
st = os.stat(sh_file_path)
os.chmod(sh_file_path, st.st_mode | stat.S_IEXEC)


sh_content = f"""#!/bin/bash
# Generated script by 1_create_sh_script.py
DATABASES_DIR="inputs/databases/"
REF_SPECIES_PATH="inputs/reference/"
REF_SPECIES_NO_EXT_PATH="inputs/reference/"
NONREF_SPECIES_DIR="{config['query_proteins_dir']}"
OUTPUT_DIR="outputs/"

# Find the file in the directory that does not have the .dmnd extension
for file in "$REF_SPECIES_PATH"/*; do
    if [[ -f "$file" && "$file" != *.dmnd ]]; then
        # Extract the filename
        filename_with_ext=$(basename "$file")
        REF_SPECIES_PATH=$REF_SPECIES_PATH$filename_with_ext
        REF_SPECIES_NO_EXT_PATH=$REF_SPECIES_NO_EXT_PATH${{filename_with_ext%.*}}
        break
    fi
done

echo "Query Non-Ref vs Ref"
for FILE in $(ls $NONREF_SPECIES_DIR)
do
    SPECIES_NAME="$(echo "$FILE" | cut -d'.' -f1)"
    OUT=${{OUTPUT_DIR}}${{SPECIES_NAME}}_vs_ref.tsv

    diamond blastp -q $NONREF_SPECIES_DIR$FILE -d $REF_SPECIES_NO_EXT_PATH -o $OUT --ultra-sensitive --top 0 --quiet
done

echo "Query Ref vs Non-Ref"
for DB in $(ls $DATABASES_DIR)
do
    SPECIES_NAME="$(echo "$DB" | cut -d'.' -f1)"
    DB_NAME="${{DATABASES_DIR}}${{SPECIES_NAME}}"
    OUT="${{OUTPUT_DIR}}ref_vs_${{SPECIES_NAME}}.tsv"

    diamond blastp -q $REF_SPECIES_PATH -d $DB_NAME -o $OUT --ultra-sensitive --top 0 --quiet
done
echo "Done."
"""

sh_file_path = '5_run_bidirectional_diamond.sh'

with open(sh_file_path, 'w') as file:
    file.write(sh_content)

# Give exec permissions
st = os.stat(sh_file_path)
os.chmod(sh_file_path, st.st_mode | stat.S_IEXEC)