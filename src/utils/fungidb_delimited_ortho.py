import re
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def fungidb_delimited_ortho(row):
    ortho_f_name = row['cds_file_name'].replace('_from_gff', '')

    ortho_genes = list()
    try:
        ortho_genes = list(SeqIO.parse('src/utils/ortholog_finder/orthogroups/' + ortho_f_name, 'fasta'))
    except:
        return row['original_name']

    if len(ortho_genes) == 0:
        return row['original_name']
    
    pids = list()
    descriptions = list()
    ortho_to_gene = list()
    ortho_to_ref_prot = list()
    ref_species = list()

    for o_gene in ortho_genes:
        o_gene_pid_match = re.search(r'\[protein_id=(.*?)\]', o_gene.description)
        ortho_to_gene_match = re.search(r'\[orthologous_to_gene=(.*?)\]', o_gene.description)
        ortho_to_ref_prot_match = re.search(r'\[orthologous_to_ref_protein=(.*?)\]', o_gene.description)
        ref_species_match = re.search(r'\[ref_species=(.*?)\]', o_gene.description)

        if o_gene_pid_match and ortho_to_gene_match and ortho_to_ref_prot_match and ref_species_match:
            pids.append(o_gene_pid_match.group(1))
            ortho_to_gene.append(ortho_to_gene_match.group(1))
            ortho_to_ref_prot.append(ortho_to_ref_prot_match.group(1))
            ref_species.append(ref_species_match.group(1))
            descriptions.append(o_gene.description)
        else:
            print('Problem in ', ortho_f_name, o_gene)
            return row['original_name']

    genes = SeqIO.parse('data/fourdbs_concat/delimited_cds_from_gff/' + row['cds_file_name'], 'fasta')

    count = 0
    found = list()
    with open('data/fourdbs_concat/ortho_from_gff/' + ortho_f_name, 'w') as file:
        for g in list(genes):
            gene_pid_match = g.id
            
            if gene_pid_match in pids:
                count += 1
                found.append(gene_pid_match)

                idx = pids.index(gene_pid_match)
                rec = SeqRecord(
                    g.seq,
                    id=gene_pid_match,
                    description=(descriptions[idx]),
                )
                SeqIO.write(rec, file, 'fasta')

    if count != len(pids):
        print(f'Problem in {ortho_f_name} - Found {count} genes in delimited_cds_from_gff but {len(pids)} genes in orthogroups')
        print(f'{found} =/= {pids}')
        return row['original_name']
    
    return 'OK'
            