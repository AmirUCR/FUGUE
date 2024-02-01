import os
import shutil
import pandas as pd


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
            'mucor_lanceolatus', 'mucor_fuscus', 'mucor_endophyticus']

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

    if not os.path.exists(new_cds_dir): os.makedirs(new_cds_dir)

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

    for idx, row in concat.iterrows():
        if idx % 100 == 0:
            print(f'Done with {idx} species...', end='\r')

        cds_f_name = row['cds_file_name']

        shutil.copy(f"{source_to_dir[row['source']]}/cds_from_gff/{cds_f_name}", new_cds_dir)

    concat.to_csv(os.path.join(concat_destination_dir, 'fourdbs_gff_input_species.csv'), index=False)
    print()
    print(f"Done. Check {os.path.join(concat_destination_dir, 'fourdbs_gff_input_species.csv')}")