import os
import pandas as pd

acc_to_name = {
    'NP_009673.1': 'LYS2',
    'NP_010290.3': 'TRP1',
    'NP_010851.1': 'CAN1',
    'NP_010893.3': 'URA3',
    'NP_012965.3': 'GAP1',
    'NP_013406.1': 'MET17',
    'NP_015387.1': 'FCY1',
}

num_files = len(os.listdir('outputs'))

ortho = pd.read_csv("Orthogroups.tsv", sep='\t').set_index('Orthogroup')

def find_cols_with_x_not_found(df, x):
    return df.columns[df.eq('Not found').sum() == x].tolist()

# Identify columns where only n row is "Not Found"
cols_with_one_not_found = find_cols_with_x_not_found(ortho, 1)
cols_with_two_not_found = find_cols_with_x_not_found(ortho, 2)
cols_with_three_not_found = find_cols_with_x_not_found(ortho, 3)
cols_with_four_not_found = find_cols_with_x_not_found(ortho, 4)
cols_with_five_not_found = find_cols_with_x_not_found(ortho, 5)
cols_with_six_not_found = find_cols_with_x_not_found(ortho, 6)
cols_with_seven_not_found = find_cols_with_x_not_found(ortho, 7)


d = {
    'NP_009673.1': [],  # LYS2
    'NP_010290.3': [],  # TRP1
    'NP_010851.1': [],  # CAN1
    'NP_010893.3': [],  # URA3
    'NP_012965.3': [],  # GAP1
    'NP_013406.1': [],  # MET17
    'NP_015387.1': []   # FCY1
}

missing_in = {
    'NP_009673.1': [],  # LYS2
    'NP_010290.3': [],  # TRP1
    'NP_010851.1': [],  # CAN1
    'NP_010893.3': [],  # URA3
    'NP_012965.3': [],  # GAP1
    'NP_013406.1': [],  # MET17
    'NP_015387.1': []   # FCY1
}

for species in [e for e in ortho.columns.tolist()]:
    df = pd.read_csv(f'outputs/{species}.tsv', sep='\t', names=['scer_gene', 'ortho', 'identity', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'])

    for gene, l in d.items():
        if ortho[species][gene] != 'Not found' and gene in df['scer_gene'].values:
            d[gene].append(df[df['scer_gene'] == gene]['identity'].values[0])
        else:
            missing_in[gene].append(species)
            
total_missing = 0
for gene, missing_in_species in missing_in.items():
    length = len(missing_in_species)
    print(f'Gene {gene} {acc_to_name[gene]} is missing its ortholog in {length} species')
    total_missing += length

print(f'Total missing orthologs {total_missing}')


for gene in d.keys():
    low_quality = len([e for e in d[gene] if e >= 0 and e < 30])
    length = len(d[gene])

    print(f'Gene {gene} has {low_quality}/{length} low quality (< 30% identity) orthologs')


ortho_high_identity = ortho.copy(deep=True)

for species in ortho.columns.tolist():
    df = pd.read_csv(f'outputs/{species}.tsv', sep='\t', names=['scer_gene', 'ortho', 'identity', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'])
    low_identities = df[df['identity'] < 30]['scer_gene'].tolist()
    ortho_high_identity.loc[ortho_high_identity.index.isin(low_identities), species] = 'Not found'

cols_with_six_not_found = find_cols_with_x_not_found(ortho_high_identity, 6)
cols_with_seven_not_found = find_cols_with_x_not_found(ortho_high_identity, 7)

ortho_high_identity.drop(cols_with_six_not_found, axis=1, inplace=True)
ortho_high_identity.drop(cols_with_seven_not_found, axis=1, inplace=True)

ortho_high_identity = ortho_high_identity.reset_index()
ortho_high_identity.to_csv('Orthogroups_HI.tsv', sep='\t', index=False)

ortho_high_identity = pd.read_csv('Orthogroups_HI.tsv', sep='\t').set_index('Orthogroup')

d = {
    'NP_009673.1': [],  # LYS2
    'NP_010290.3': [],  # TRP1
    'NP_010851.1': [],  # CAN1
    'NP_010893.3': [],  # URA3
    'NP_012965.3': [],  # GAP1
    'NP_013406.1': [],  # MET17
    'NP_015387.1': []   # FCY1
}

missing_in = {
    'NP_009673.1': [],  # LYS2
    'NP_010290.3': [],  # TRP1
    'NP_010851.1': [],  # CAN1
    'NP_010893.3': [],  # URA3
    'NP_012965.3': [],  # GAP1
    'NP_013406.1': [],  # MET17
    'NP_015387.1': []   # FCY1
}

for species in [e for e in ortho_high_identity.columns.tolist()]:
    df = pd.read_csv(f'outputs/{species}.tsv', sep='\t', names=['scer_gene', 'ortho', 'identity', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'])

    for gene, l in d.items():
        if ortho_high_identity[species][gene] != 'Not found' and gene in df['scer_gene'].values:
            d[gene].append(df[df['scer_gene'] == gene]['identity'].values[0])
        else:
            missing_in[gene].append(species)
            
total_missing = 0
for gene, missing_in_species in missing_in.items():
    length = len(missing_in_species)
    print(f'Gene {gene} {acc_to_name[gene]} is missing its ortholog in {length} species')
    total_missing += length

print(f'Total missing orthologs {total_missing}')

ortho_high_identity = ortho_high_identity.reset_index()

import os

if not os.path.exists('../../../data/fourdbs_concat/removed'):
    os.makedirs('../../../data/fourdbs_concat/removed/cds')
    os.makedirs('../../../data/fourdbs_concat/removed/genomes')
    os.makedirs('../../../data/fourdbs_concat/removed/proteomes')
    os.makedirs('../../../data/fourdbs_concat/removed/cds_from_gff')

input_species = pd.read_csv('../../../data/fourdbs_concat/fourdbs_input_species.csv')
ortho_species = ortho_high_identity.columns.tolist()

removed_species = input_species[~input_species['original_name'].isin(ortho_species)].reset_index(drop=True)
removed_species.to_csv('../../../data/fourdbs_concat/removed/removed_species.csv', index=False)

input_species = input_species[input_species['original_name'].isin(ortho_species)].reset_index(drop=True)
input_species.to_csv('../../../data/fourdbs_concat/fourdbs_input_species.csv', index=False)

concat_destination_dir = '../../../data/fourdbs_concat'
cds = os.listdir(concat_destination_dir + "/cds/")
gffs = os.listdir(concat_destination_dir + "/gff/")
genomes = os.listdir(concat_destination_dir + "/genomes/")
proteomes = os.listdir(concat_destination_dir + "/proteomes/")
cds_from_gffs = os.listdir(concat_destination_dir + "/cds_from_gff/")
delim_cds_from_gffs = os.listdir(concat_destination_dir + "/delimited_cds_from_gff/")

for cds_f in cds:
    if cds_f.replace('_cds', '_cds_from_gff') not in input_species['cds_file_name'].tolist():
        os.remove(f'{concat_destination_dir}/cds/{cds_f}')
        print('remove', f'{concat_destination_dir}/cds/{cds_f}')

for gff_f in gffs:
    if gff_f not in input_species['gff_file_name'].tolist():
        os.remove(f'{concat_destination_dir}/gff/{gff_f}')
        print('remove', f'{concat_destination_dir}/gff/{gff_f}')

for genome_f in genomes:
    if genome_f not in input_species['genome_file_name'].tolist():
        os.remove(f'{concat_destination_dir}/genomes/{genome_f}')
        print('remove', f'{concat_destination_dir}/genomes/{genome_f}')

for prot_f in proteomes:
    if prot_f.replace('.faa', '') not in input_species['original_name'].tolist():
        os.remove(f'{concat_destination_dir}/proteomes/{prot_f}')
        print('remove', f'{concat_destination_dir}/proteomes/{prot_f}')

for cds_from_gff in cds_from_gffs:
    if cds_from_gff not in input_species['cds_file_name'].tolist():
        os.remove(f'{concat_destination_dir}/cds_from_gff/{cds_from_gff}')
        print('remove', f'{concat_destination_dir}/cds_from_gff/{cds_from_gff}')

for delim_cds_from_gff in delim_cds_from_gffs:
    if delim_cds_from_gff not in input_species['cds_file_name'].tolist():
        os.remove(f'{concat_destination_dir}/delimited_cds_from_gff/{delim_cds_from_gff}')
        print('remove', f'{concat_destination_dir}/delimited_cds_from_gff/{delim_cds_from_gff}')