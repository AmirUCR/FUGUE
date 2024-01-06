import re
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def fungidb_gff_to_cds(name, cds_file, genome_file, gff_file, output_file):
    # load the CDS we are interested in aka orthologs
    records = SeqIO.parse(cds_file, 'fasta')

    goi = list()
    headers = dict()
    for record in records:
        gene_id = record.id
        goi.append(gene_id)
        headers[gene_id] = record.description

    # Load the genome file
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

    # Function to extract CDS from GFF3
    def extract_cds(gff_file, genome, goi_list):
        cds_sequences = {key: '' for key in goi_list}
        reading_gene = ''
        current_seq_id = ''
        current_gene_strand = ''
        current_coords: list[tuple[str, str]] = list()

        with open(gff_file, 'r') as gff:
            for line in gff:
                if line.startswith('#'):
                    reading_gene = ''
                    current_seq_id = ''
                    current_gene_strand = ''
                    current_coords = []
                    continue
                parts = line.strip().split('\t')

                if len(parts) < 2:
                    print('fungidb_gff_to_csv.py: Problematic line. Skipping')
                    print(name, line)
                    continue

                if parts[2].lower() == 'protein_coding_gene':
                    # process previous gene
                    current_coords = sorted(current_coords, key=lambda x: x[0])

                    for start, end in current_coords:
                        seq = genome[current_seq_id].seq[start-1:end]
                        if current_gene_strand == '-':
                            seq = seq.reverse_complement()
                            cds_sequences[reading_gene] = str(seq) + '|' + cds_sequences[reading_gene]
                        else:
                            cds_sequences[reading_gene] += str(seq) + '|'

                    reading_gene = ''
                    current_seq_id = ''
                    current_gene_strand = ''
                    current_coords = []

                    gene_id = re.search('ID=([^;]+)', parts[8])
                    
                    if gene_id:
                        gene_id = gene_id.group(1)
                    else:
                        print(f'fungidb_gff_to_cds: {name} has no ID= in its GFF. Skipping.')
                        return
                    
                    for g in goi:
                        if g.startswith(gene_id):
                            reading_gene = g

                elif reading_gene != '' and parts[2].lower() == 'cds':
                    seq_id, start, end, strand = parts[0], int(parts[3]), int(parts[4]), parts[6]

                    current_seq_id = seq_id
                    current_gene_strand = strand
                    current_coords.append((start, end))

        return cds_sequences

    goi_seq_dict = extract_cds(gff_file, genome, goi)

    if goi_seq_dict is None:
        return name

    records = list()
    for k, v in goi_seq_dict.items():
        if v != '':
            header = headers[k]
            record = SeqRecord(
                Seq(v),
                id=header.split(' ')[0],
                description=' '.join(headers[k].split(' ')[1:]),
            )
            records.append(record)

    # Failure
    if len(records) == 0:
        return name

    new_name = name + '_cds_from_gff.fna'
    out = os.path.join(output_file, new_name)

    with open(out, "w") as f:
        SeqIO.write(records, f, "fasta")

    return 'OK'