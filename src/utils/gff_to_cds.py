from ncbi_gff_to_cds import ncbi_gff_to_cds
from fungidb_gff_to_cds import fungidb_gff_to_cds
from ensembl_gff_to_cds import ensembl_gff_to_cds
from mycocosm_gff_to_cds import mycocosm_gff_to_cds

import os
import pandas as pd

def create_cds_from_gff():
    if not os.path.exists('../../data/fourdbs_concat/cds_from_gff'):
        os.mkdir('../../data/fourdbs_concat/cds_from_gff')

    out = '../../data/fourdbs_concat/cds_from_gff'

    df = pd.read_csv('../../data/fourdbs_concat/fourdbs_hi_input_species.csv')

    ncbi_base_path = '../../data/NCBI'
    fdb_base_path = '../../data/FungiDB'
    ensembl_base_path = '../../data/EnsemblFungi'
    mycocosm_base_path = '../../data/MycoCosm'

    losers = list()

    for idx, row in df.iterrows():
        src = row['source']
        species_name = row['species_name']
        cds_file_name = row['cds_file_name']
        genome_file_name = row['genome_file_name']
        gff_file_name = row['gff_file_name']

        if os.path.exists(os.path.join(out, species_name + '_cds_from_gff.fna')):
            print(os.path.join(out, species_name + '_cds_from_gff.fna'), 'exists. Skipping.')
            continue

        response = ''

        if src == 'NCBI':
            response = ncbi_gff_to_cds(species_name,
                                    f'ortholog_finder/orthogroups/{cds_file_name}',
                                    f'{ncbi_base_path}/genomes/{genome_file_name}',
                                    f'{ncbi_base_path}/gff/{gff_file_name}',
                                    out)
        elif src == 'FungiDB':
            response = fungidb_gff_to_cds(species_name,
                                        f'ortholog_finder/orthogroups/{cds_file_name}',
                                        f'{fdb_base_path}/genomes/{genome_file_name}',
                                        f'{fdb_base_path}/gff/{gff_file_name}',
                                        out)
        elif src == 'EnsemblFungi':
            response = ensembl_gff_to_cds(species_name,
                                        f'ortholog_finder/orthogroups/{cds_file_name}',
                                        f'{ensembl_base_path}/genomes/{genome_file_name}',
                                        f'{ensembl_base_path}/gff/{gff_file_name}',
                                        out)
        elif src == 'MycoCosm':
            response = mycocosm_gff_to_cds(species_name,
                                        f'ortholog_finder/orthogroups/{cds_file_name}',
                                        f'{mycocosm_base_path}/genomes/{genome_file_name}',
                                        f'{mycocosm_base_path}/gff/{gff_file_name}',
                                        out)
        
        if response != 'OK':
            losers.append(species_name)

    print('Losers:', losers)

    df = df[~df['species_name'].isin(losers)]
    df['cds_file_name'] = df['cds_file_name'].str.replace('_cds', '_cds_from_gff')
    df.to_csv('../../data/fourdbs_concat/fourdbs_hi_gff_input_species.csv', index=False)

    print('Created CDS files under data/fourdbs_concat/cds_from_gff. Place these in the input folder of ALLEGRO and point to the folder in config.yaml')
    print('Created data/fourdbs_concat/fourdbs_hi_gff_input_species.csv Use this as the input species csv for ALLEGRO.')