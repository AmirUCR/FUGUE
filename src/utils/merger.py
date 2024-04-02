import os
import shutil
import pandas as pd
from Bio import SeqIO


def merge_dbs():
    ncbi_dir = 'data/NCBI'
    fdb_dir = 'data/FungiDB'
    mycocosm_dir = 'data/MycoCosm'
    ensembl_dir = 'data/EnsemblFungi'
    concat_destination_dir = 'data/fourdbs_concat'

    new_cds_dir = os.path.join(concat_destination_dir, 'cds')
    new_genome_dir = os.path.join(concat_destination_dir, 'genomes')
    new_prot_dir = os.path.join(concat_destination_dir, 'proteomes')
    new_gff_dir = os.path.join(concat_destination_dir, 'gff')

    if not os.path.exists(concat_destination_dir): os.mkdir(concat_destination_dir)
    if not os.path.exists(os.path.join(concat_destination_dir, 'cds')): os.makedirs(new_cds_dir)
    if not os.path.exists(os.path.join(concat_destination_dir, 'genomes')): os.makedirs(new_genome_dir)
    if not os.path.exists(os.path.join(concat_destination_dir, 'proteomes')): os.makedirs(new_prot_dir)
    if not os.path.exists(os.path.join(concat_destination_dir, 'gff')): os.makedirs(new_gff_dir)

    ncbi = pd.read_csv(f'{ncbi_dir}/ncbi_input_species.csv')
    ncbi['source'] = 'NCBI'

    fdb = pd.read_csv(f'{fdb_dir}/fungidb_input_species.csv')
    fdb['source'] = 'FungiDB'

    mycocosm = pd.read_csv(f'{mycocosm_dir}/mycocosm_input_species.csv')
    mycocosm['source'] = 'MycoCosm'

    ensembl = pd.read_csv(f'{ensembl_dir}/ensemblfungi_input_species.csv')
    ensembl['source'] = 'EnsemblFungi'

    source_to_dir = {
        'NCBI': ncbi_dir,
        'FungiDB': fdb_dir,
        'MycoCosm': mycocosm_dir,
        'EnsemblFungi': ensembl_dir
    }

    fdb = fdb.drop_duplicates(subset='species_name').reset_index(drop=True)
    ncbi = ncbi.drop_duplicates(subset='species_name').reset_index(drop=True)
    ensembl = ensembl.drop_duplicates(subset='species_name').reset_index(drop=True)
    mycocosm = mycocosm.drop_duplicates(subset='species_name').reset_index(drop=True)

    concat = pd.concat([ncbi, fdb, ensembl, mycocosm], ignore_index=True)
    concat = concat.drop_duplicates(subset=['species_name'])

    # Losers go here
    losers = ['candida_glabrata', 'lasallia_pustulata', 'neosartorya_fischeri_nrrl_181',
            'aspergillus_campestris', 'volvariella_volvacea_v23', 'mucor_racemosus',
            'mucor_lanceolatus', 'mucor_fuscus', 'mucor_endophyticus',
            'puccinia_coronata_var_avenae_f_sp_avenae']

    concat = concat[~concat['original_name'].isin(losers)]
    concat = concat.reset_index(drop=True)

    for idx, row in concat.iterrows():
        if idx % 100 == 0: print(idx, '...')

        cds_f_name = row['cds_file_name']
        genome_f_name = row['genome_file_name']
        protein_f_name = row['original_name'] + '.faa'
        gff_f_name = row['gff_file_name']

        shutil.copy(f"{source_to_dir[row['source']]}/cds/{cds_f_name}", new_cds_dir)
        shutil.copy(f"{source_to_dir[row['source']]}/genomes/{genome_f_name}", new_genome_dir)
        shutil.copy(f"{source_to_dir[row['source']]}/proteomes/{protein_f_name}", new_prot_dir)
        shutil.copy(f"{source_to_dir[row['source']]}/gff/{gff_f_name}", new_gff_dir)

    concat.to_csv(os.path.join(concat_destination_dir, 'fourdbs_input_species.csv'), index=False)


def merge_gffs():
    ncbi_dir = 'data/NCBI'
    fdb_dir = 'data/FungiDB'
    mycocosm_dir = 'data/MycoCosm'
    ensembl_dir = 'data/EnsemblFungi'
    concat_destination_dir = 'data/fourdbs_concat'

    new_cds_dir = os.path.join(concat_destination_dir, 'cds_from_gff')
    new_delimited_cds_dir = os.path.join(concat_destination_dir, 'delimited_cds_from_gff')

    if not os.path.exists(new_cds_dir):
        os.makedirs(new_cds_dir)

    if not os.path.exists(new_delimited_cds_dir):
        os.makedirs(new_delimited_cds_dir)

    ncbi = pd.read_csv(f'{ncbi_dir}/ncbi_gff_input_species.csv')
    ncbi['source'] = 'NCBI'

    fdb = pd.read_csv(f'{fdb_dir}/fungidb_gff_input_species.csv')
    fdb['source'] = 'FungiDB'

    mycocosm = pd.read_csv(f'{mycocosm_dir}/mycocosm_gff_input_species.csv')
    mycocosm['source'] = 'MycoCosm'

    ensembl = pd.read_csv(f'{ensembl_dir}/ensemblfungi_gff_input_species.csv')
    ensembl['source'] = 'EnsemblFungi'

    source_to_dir = {
        'NCBI': ncbi_dir,
        'FungiDB': fdb_dir,
        'MycoCosm': mycocosm_dir,
        'EnsemblFungi': ensembl_dir
    }

    fdb = fdb.drop_duplicates(subset='species_name').reset_index(drop=True)
    ncbi = ncbi.drop_duplicates(subset='species_name').reset_index(drop=True)
    ensembl = ensembl.drop_duplicates(subset='species_name').reset_index(drop=True)
    mycocosm = mycocosm.drop_duplicates(subset='species_name').reset_index(drop=True)

    concat = pd.concat([ncbi, fdb, ensembl, mycocosm], ignore_index=True)
    concat = concat.drop_duplicates(subset=['species_name'])

    # Losers go here
    losers = ['candida_glabrata', 'lasallia_pustulata', 'neosartorya_fischeri_nrrl_181',
            'aspergillus_campestris', 'volvariella_volvacea_v23', 'mucor_racemosus',
            'mucor_lanceolatus', 'mucor_fuscus', 'mucor_endophyticus']

    concat = concat[~concat['original_name'].isin(losers)]
    concat = concat.reset_index(drop=True)

    losers = list()

    for idx, row in concat.iterrows():
        if idx % 100 == 0:
            print(f'Done with {idx} species...', end='\r')

        cds_f_name = row['cds_file_name']
        
        l = len(list(SeqIO.parse(f"{source_to_dir[row['source']]}/cds_from_gff/{cds_f_name}", 'fasta')))

        if l < 2000:
            losers.append(row['original_name'])
            print(f'merger.py: File {cds_f_name} has fewer than 2000 genes ({l}). Skipping.')
            continue

        shutil.copy(f"{source_to_dir[row['source']]}/cds_from_gff/{cds_f_name}", new_cds_dir)

        if os.path.exists(f"{source_to_dir[row['source']]}/delimited_cds_from_gff/{cds_f_name}"):
            shutil.copy(f"{source_to_dir[row['source']]}/delimited_cds_from_gff/{cds_f_name}", new_delimited_cds_dir)

    concat = concat[~concat['original_name'].isin(losers)]
    concat = concat.reset_index(drop=True)

    cds = os.listdir(concat_destination_dir + "/cds/")
    gffs = os.listdir(concat_destination_dir + "/gff/")
    genomes = os.listdir(concat_destination_dir + "/genomes/")
    proteomes = os.listdir(concat_destination_dir + "/proteomes/")
    cds_from_gffs = os.listdir(concat_destination_dir + "/cds_from_gff/")

    remove_list = list()
    for cds_f in concat['cds_file_name'].tolist():
        if cds_f.replace('_cds_from_gff', '_cds') not in cds:
            remove_list.append(cds_f)
    concat = concat[~concat['cds_file_name'].isin(remove_list)].reset_index(drop=True)

    for cds_f in cds:
        if cds_f.replace('_cds.fna', '_cds_from_gff.fna') not in concat['cds_file_name'].tolist():
            os.remove(f'{concat_destination_dir}/cds/{cds_f}')
            print('remove', f'{concat_destination_dir}/cds/{cds_f}')

    for gff_f in gffs:
        if gff_f not in concat['gff_file_name'].tolist():
            os.remove(f'{concat_destination_dir}/gff/{gff_f}')
            print('remove', f'{concat_destination_dir}/gff/{gff_f}')

    for genome_f in genomes:
        if genome_f not in concat['genome_file_name'].tolist():
            os.remove(f'{concat_destination_dir}/genomes/{genome_f}')
            print('remove', f'{concat_destination_dir}/genomes/{genome_f}')

    for prot_f in proteomes:
        if prot_f.replace('.faa', '') not in concat['original_name'].tolist():
            os.remove(f'{concat_destination_dir}/proteomes/{prot_f}')
            print('remove', f'{concat_destination_dir}/proteomes/{prot_f}')

    for cds_from_gff in cds_from_gffs:
        if cds_from_gff not in concat['cds_file_name'].tolist():
            os.remove(f'{concat_destination_dir}/cds_from_gff/{cds_from_gff}')
            print('remove', f'{concat_destination_dir}/cds_from_gff/{cds_from_gff}')

    concat.to_csv(os.path.join(concat_destination_dir, 'fourdbs_input_species.csv'), index=False)

    print()
    print(f'Merged {len(concat)} species to {concat_destination_dir}/fourdbs_input_species/.')
    print(f"Done. Check {os.path.join(concat_destination_dir, 'fourdbs_input_species.csv')}")