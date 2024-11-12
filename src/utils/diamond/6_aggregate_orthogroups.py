import os
import pandas as pd
from Bio import SeqIO 

outputs_directory = 'outputs'

headers = ['query_accession', 'target_accession', 'sequence_identity', 'length', 'mismatches',
           'gap_openings', 'query_start', 'query_end', 'target_start', 'target_end', 'e_value',
           'bit_score']

ref_species = ''
for f in os.listdir('inputs/reference'):
    if not f.endswith('.dmnd'):
        ref_species = f

ref_records = SeqIO.parse(f'inputs/reference/{ref_species}', 'fasta')
ref_gene_accessions = list()

for rec in ref_records:
    ref_gene_accessions.append(rec.id)

orthogroups_dict = {
    'Orthogroup': ref_gene_accessions
    }

empty_files = list()
file_pairs = list()

for f in os.listdir('outputs'):
    if os.stat("outputs/" + f ).st_size == 0:
        if f.endswith('_vs_ref.tsv'):
            empty_files.append(f)
            empty_files.append('ref_vs_' + f.split('_vs_ref.tsv')[0] + '.tsv')

for f in os.listdir('outputs'):
    if f.endswith('_vs_ref.tsv') and f not in empty_files:
        file_pairs.append((f, 'ref_vs_' + f.split('_vs_ref.tsv')[0] + '.tsv'))

def load_blast_results(file):
    return pd.read_csv(file, sep='\t', header=None,
                       names=['query', 'subject', 'identity', 'length', 'mismatches', 'gap_openings', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score'])

# Function to get the best hit based on the highest bit score (or lowest e-value)
def get_best_hits(df):
    return df.loc[df.groupby('query')['bit_score'].idxmax()]

for pair in file_pairs:
    nrs_vs_rs_tsv = pair[0]  # 'sclerotinia_borealis_vs_ref.tsv'
    rs_vs_nrs_tsv = pair[1]  # 'ref_vs_sclerotinia_borealis.tsv'

    nrs_name = nrs_vs_rs_tsv.split('_vs_ref.tsv')[0]
    orthogroups_dict[nrs_name] = list()

    nrs_vs_rs_path = os.path.join(outputs_directory, nrs_vs_rs_tsv)
    rs_vs_nrs_path = os.path.join(outputs_directory, rs_vs_nrs_tsv)

    other_vs_ref = load_blast_results(nrs_vs_rs_path)
    ref_vs_other = load_blast_results(rs_vs_nrs_path)

    # Get best hits
    sac_vs_other_best = get_best_hits(ref_vs_other)
    other_vs_sac_best = get_best_hits(other_vs_ref)

    # Merge results with suffixes to differentiate columns
    merged_results = sac_vs_other_best.merge(other_vs_sac_best, left_on=['query', 'subject'], right_on=['subject', 'query'], suffixes=('_sac', '_other'))

    # Extract reciprocal best hits
    nrs_df = merged_results[['query_sac', 'subject_sac', 'identity_sac', 'length_sac', 'mismatches_sac', 'gap_openings_sac', 'q_start_sac', 'q_end_sac', 's_start_sac', 's_end_sac', 'evalue_sac', 'bit_score_sac']]
    nrs_df.columns = ['query', 'subject', 'identity', 'length', 'mismatches', 'gap_openings', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score']

    nrs_df.to_csv(f'outputs/reciprocal/{nrs_name}.tsv', index=False, sep='\t')

    for ref_gene_accession in ref_gene_accessions:
        target_accession = 'Not found'
        if ref_gene_accession in nrs_df['query'].values:
            target_accession = nrs_df[nrs_df['query'] == ref_gene_accession]['subject'].tolist()[0]

            if len(nrs_df[nrs_df['query'] == ref_gene_accession]['subject'].tolist()) > 1:
                print(nrs_name)

        orthogroups_dict[nrs_name].append(target_accession)

pd.DataFrame.from_dict(orthogroups_dict).to_csv('Orthogroups.tsv', sep='\t', index=False)