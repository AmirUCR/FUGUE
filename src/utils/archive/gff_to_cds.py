from utils.ncbi_gff_to_cds import ncbi_gff_to_cds
from utils.fungidb_gff_to_cds import fungidb_gff_to_cds
from utils.ensembl_gff_to_cds import ensembl_gff_to_cds
from utils.mycocosm_gff_to_cds import mycocosm_gff_to_cds

import os
import pandas as pd


# ALL PATHS ARE RELATIVE TO THE ROOT DIRECTORY
def create_cds_from_gff():
    sources = ['NCBI', 'FungiDB', 'EnsemblFungi', 'MycoCosm']
    base_paths = ['data/NCBI', 'data/FungiDB', 'data/EnsemblFungi', 'data/MycoCosm']
    manifests = ['ncbi_input_species.csv', 'fungidb_input_species.csv', 'ensemblfungi_input_species.csv', 'mycocosm_input_species.csv']

    losers = list()

    count = 0
    for m_idx, m in enumerate(manifests):
        csv_path = os.path.join(base_paths[m_idx], m)

        df = pd.read_csv(csv_path)
        for _, row in df.iterrows():
            if count % 100 == 0:
                print(f'Done with {count} species...')

            src = sources[m_idx]
            src_base_path = base_paths[m_idx]

            original_name = row['original_name']
            cds_file_name = row['cds_file_name']
            genome_file_name = row['genome_file_name']
            gff_file_name = row['gff_file_name']

            if os.path.exists(os.path.join(f'data/{src}/cds_from_gff/', original_name + '_cds_from_gff.fna')):
                print(os.path.join(f'data/{src}/cds_from_gff/', original_name + '_cds_from_gff.fna'), 'exists. Skipping.')
                count += 1
                continue

            response = ''

            if src == 'NCBI':
                response = ncbi_gff_to_cds(original_name,
                                        f'{src_base_path}/cds/{cds_file_name}',
                                        f'{src_base_path}/genomes/{genome_file_name}',
                                        f'{src_base_path}/gff/{gff_file_name}')
            elif src == 'FungiDB':
                response = fungidb_gff_to_cds(original_name,
                                            f'{src_base_path}/cds/{cds_file_name}',
                                            f'{src_base_path}/genomes/{genome_file_name}',
                                            f'{src_base_path}/gff/{gff_file_name}')
            elif src == 'EnsemblFungi':
                response = ensembl_gff_to_cds(original_name,
                                            f'{src_base_path}/cds/{cds_file_name}',
                                            f'{src_base_path}/genomes/{genome_file_name}',
                                            f'{src_base_path}/gff/{gff_file_name}')
            elif src == 'MycoCosm':
                response = mycocosm_gff_to_cds(original_name,
                                            f'{src_base_path}/cds/{cds_file_name}',
                                            f'{src_base_path}/genomes/{genome_file_name}',
                                            f'{src_base_path}/gff/{gff_file_name}')
            
            if response != 'OK':
                print(f'{src}: No records found for {src_base_path}/cds/{cds_file_name}. Skipping.')
                losers.append(original_name)
            
            count += 1
        
        df = df[~df['original_name'].isin(losers)]
        df['cds_file_name'] = df['cds_file_name'].str.replace('_cds', '_cds_from_gff')
        df.to_csv(f'data/{src}/{src.lower()}_gff_input_species.csv', index=False)

    print()
    print('Losers:', losers)
