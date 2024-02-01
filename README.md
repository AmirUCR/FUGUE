# FUGUE / Fungal Universe Genome Unification Engine

This repository contains instructions and scripts that allow you to download genomes, CDS, GFF, and proteomes for more than 2,000 fungal species used in ALLEGRO.

## Prerequisites
1. First, download Miniconda [https://docs.conda.io/en/main/miniconda.html](https://docs.conda.io/en/main/miniconda.html)
2. Clone this repository by either clicking on the green Code on the top right and clicking "Download ZIP," or downloading [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) and then `$ git clone https://github.com/AmirUCR/Fugue.git` in your desired directory.
3. These installation instructions are for Ubuntu 20.04.6 LTS. Create a conda environment and activate it:

```
conda create -n fugue python=3.10 -y
conda activate fugue
```

4. Install the required Python packages.

```
pip install pandas
pip install biopython
pip install PyYAML
pip install requests
```

## Gathering the Data
The order you should download from these databases is as follows:

1. NCBI
2. FungiDB
3. EnsemblFungi
4. MycoCosm

For steps 1 through 4, visit here: https://drive.google.com/drive/folders/1FSkpgUBtfJ4NcyftYYQicKJLKNVGUR9d?usp=drive_link


You need to perform the mandatory steps in [Data Source] NCBI Datasets and [Data Source] MycoCosm. These two databases need special authentication methods AKA your own credentials. The other two documents for data sources are for your reference and no further action is needed from you.


Run `$ python src/main.py` and download 1 through 4. Then merge the databases using option 5. 

In ALLEGRO, we run DIAMOND to find orthogroups of *S. cerevisiae* genes of interest across all other species. Navigate to src/utils/diamond. Ensure you have the diamond python package installed in a conda environment called 'diamond' (alternatively, modify the name of the environment in the shell script 0_just_run_this.sh). Ensure you have the allegro environment installed (it must have biopython, pyyaml, and pandas installed). Place the .faa amino acid file for *S. cerevisiae* (or your reference species of interest) in inputs/reference. Modify make_proteome_config.yaml to point to the correct CDS and proteome files for your reference species. Run 0_just_run_this.sh. This will generate Orthogroups.tsv. Run all cells in 5_analysis_missing_orthologs.ipynb to remove groups with < 30% protein identity and generate the new species list fourdbs_hi_input_species.csv under data/fourdbs_concat/ and generate Orthogroups_HI.tsv.

Navigate to ortholog_finder and set the parameters in the config.yaml file. Run find_orthogroup.py to generate fasta files for each species with only the genes of interest and their orthologs under orthogroups/. You may now place these files in the input directory of ALLEGRO and run the program on these ortholog genes.

If you need to use the GFF files for each species to find the intron/exon boundaries in each ortho gene, run src/main.py again and select option 6. This will generate a directory under data/fourdbs_concat/cds_from_gff/ with each orthogene having pipe | delimiters where the intron/exon boundary happens. ALLEGRO will ignore guides with pipe | in their sequences. You may use the files in cds_from_gff as input to ALLEGRO along with fourdbs_hi_gff_input_species.csv as the manifest.
