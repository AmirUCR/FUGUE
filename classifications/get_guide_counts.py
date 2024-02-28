import pandas as pd

df = pd.read_csv("species_w_classifications.csv")

subset = df[df['classification'].str.contains('Ascomycota;')].copy()
subset['species_name'] = subset['species_name'].str.replace(' ', '_').str.lower()
subset = subset.reset_index(drop=True)

input_species = pd.read_csv('../../ALLEGRO/data/input/fourdbs_gff_input_species.csv')
input_species_subset = input_species[input_species['species_name'].isin(subset['species_name'])].reset_index(drop=True)

input_species_subset.to_csv('lol.csv')

from test_guide_finder import GuideFinder
import multiprocessing
from Bio import SeqIO
import os

gf = GuideFinder()

def count_records_in_file(filename):
    print('started', filename)

    records = list(SeqIO.parse(f'{filename}', 'fasta'))
    
    num_genes = len(records)
    this_species_guides_count = 0

    for record in records:
        guides, _, _, _ = gf.identify_guides_and_indicate_strand('NGG', str(record.seq).upper(), 20, 0, 0, True)
        this_species_guides_count += len(guides)

    print('finished', filename)
    return this_species_guides_count, num_genes


def process_file(filename):
    return filename, count_records_in_file(filename)


def process_files(files):
    return [process_file(file) for file in files]


def count_records(directory: str, species_file_names: list[str]) -> int:
    # Get the list of files in the directory
    files = [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith(('.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa', '.frn')) and file in species_file_names]

    # Set the number of processes
    num_processes = len(files) // 128

    # Split the files into chunks for each process
    chunks = [files[i:i + num_processes] for i in range(0, len(files), num_processes)]

    # Create a pool of processes
    pool = multiprocessing.Pool(processes=128)

    # Map the process_files function to each chunk and get the results
    results = pool.map(process_files, chunks)

    return results


results = count_records('../../ALLEGRO/data/input/cds/cds_from_gff', input_species_subset['cds_file_name'].tolist()[:500])

import pickle

with open('results.pickle', 'wb') as f:
    pickle.dump(results, f)