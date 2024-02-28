import os
import pandas as pd
from utils.ncbi_delimited_ortho import ncbi_delimited_ortho
from utils.ensembl_delimited_ortho import ensembl_delimited_ortho
from utils.fungidb_delimited_ortho import fungidb_delimited_ortho
from utils.mycocosm_delimited_ortho import mycocosm_delimited_ortho


def delimited_ortho():
    if not os.path.exists('src/utils/ortholog_finder/orthogroups/'):
        print('Path src/utils/ortholog_finder/orthogroups/ does not exist. Exiting.')
        return
    
    if not os.path.exists('data/fourdbs_concat/ortho_from_gff'):
        os.mkdir('data/fourdbs_concat/ortho_from_gff')

    concat_destination_dir = 'data/fourdbs_concat'

    df = pd.read_csv(f'{concat_destination_dir}/fourdbs_input_species.csv')

    sources = ['NCBI', 'FungiDB', 'EnsemblFungi', 'MycoCosm']

    response = ''
    losers = list()

    for source in sources:
        slice = df[df['source'] == source]

        for _, row in slice.iterrows():
            if source == 'NCBI':
                response = ncbi_delimited_ortho(row)

            elif source == 'FungiDB':
                response = fungidb_delimited_ortho(row)

            elif source == 'EnsemblFungi':
                response = ensembl_delimited_ortho(row)

            elif source == 'MycoCosm':
                response = mycocosm_delimited_ortho(row)

            if response != 'OK':
                losers.append(response)

    print('Removed the following species:', losers)

    df = df[~df['original_name'].isin(losers)]
    df = df.reset_index(drop=True)

    cds = os.listdir(concat_destination_dir + "/cds/")
    gffs = os.listdir(concat_destination_dir + "/gff/")
    genomes = os.listdir(concat_destination_dir + "/genomes/")
    proteomes = os.listdir(concat_destination_dir + "/proteomes/")
    cds_from_gffs = os.listdir(concat_destination_dir + "/cds_from_gff/")
    delim_cds_from_gffs = os.listdir(concat_destination_dir + "/delimited_cds_from_gff/")
    ortho_from_gffs = os.listdir(concat_destination_dir + "/ortho_from_gff/")

    for cds_f in cds:
        if cds_f.replace('_cds', '_cds_from_gff') not in df['cds_file_name'].tolist():
            os.remove(f'{concat_destination_dir}/cds/{cds_f}')
            print('remove', f'{concat_destination_dir}/cds/{cds_f}')

    for gff_f in gffs:
        if gff_f not in df['gff_file_name'].tolist():
            os.remove(f'{concat_destination_dir}/gff/{gff_f}')
            print('remove', f'{concat_destination_dir}/gff/{gff_f}')

    for genome_f in genomes:
        if genome_f not in df['genome_file_name'].tolist():
            os.remove(f'{concat_destination_dir}/genomes/{genome_f}')
            print('remove', f'{concat_destination_dir}/genomes/{genome_f}')

    for prot_f in proteomes:
        if prot_f.replace('.faa', '') not in df['original_name'].tolist():
            os.remove(f'{concat_destination_dir}/proteomes/{prot_f}')
            print('remove', f'{concat_destination_dir}/proteomes/{prot_f}')

    for cds_from_gff in cds_from_gffs:
        if cds_from_gff not in df['cds_file_name'].tolist():
            os.remove(f'{concat_destination_dir}/cds_from_gff/{cds_from_gff}')
            print('remove', f'{concat_destination_dir}/cds_from_gff/{cds_from_gff}')

    for delim_cds_from_gff in delim_cds_from_gffs:
        if delim_cds_from_gff not in df['cds_file_name'].tolist():
            os.remove(f'{concat_destination_dir}/delimited_cds_from_gff/{delim_cds_from_gff}')
            print('remove', f'{concat_destination_dir}/delimited_cds_from_gff/{delim_cds_from_gff}')

    for ortho_from_gff in ortho_from_gffs:
        if ortho_from_gff.replace('_cds', '_cds_from_gff') not in df['cds_file_name'].tolist():
            os.remove(f'{concat_destination_dir}/ortho_from_gff/{ortho_from_gff}')
            print('remove', f'{concat_destination_dir}/ortho_from_gff/{ortho_from_gff}')

    df['ortho_file_name'] = df['cds_file_name'].str.replace('_cds_from_gff', '_cds')
    df = df[[c for c in df.columns if c not in ['cds_url', 'genome_url', 'proteome_url', 'gff_url', 'source']] + ['cds_url', 'genome_url', 'proteome_url', 'gff_url', 'source']]

    df.to_csv('data/fourdbs_concat/fourdbs_input_species.csv', index=False)

    print()
    print(f'Merged {len(df)} species to {concat_destination_dir}/fourdbs_input_species/.')
    print(f"Done. Check {concat_destination_dir}/fourdbs_input_species.csv")