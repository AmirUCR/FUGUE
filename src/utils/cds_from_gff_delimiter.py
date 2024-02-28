import re
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from concurrent.futures import ThreadPoolExecutor


def process_file(src, f_name):
    if 'cds_from_gff' in f_name:
        genes = SeqIO.parse(f'data/{src}/cds_from_gff/{f_name}', 'fasta')

        with open(f'data/{src}/delimited_cds_from_gff/{f_name}', 'w') as file:
            for g in genes:
                delimited_seq = ''
                
                segs = re.search(r'segs:(.*?)\s', g.description)

                if segs:
                    segs = segs.group(1)
                    for indices in segs.split(','):
                        left, right = [int(idx) for idx in indices.split('-')]
                        delimited_seq += g.seq[left - 1:right] + '|'
                else:
                    print(f'No pattern segs:(.*?)\\s found in {g.id}')
                    continue

                rec = SeqRecord(
                    seq=delimited_seq[:-1],
                    id=g.id,
                    description=(g.description),
                )
                SeqIO.write(rec, file, 'fasta')


def cds_from_gff_delimiter():
    sources = ['NCBI', 'FungiDB', 'EnsemblFungi', 'MycoCosm']

    for src in sources:
        output_dir = f'data/{src}/delimited_cds_from_gff'
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        files = [f for f in os.listdir(f'data/{src}/cds_from_gff') if 'cds_from_gff' in f]
        
        with ThreadPoolExecutor() as executor:
            futures = [executor.submit(process_file, src, f_name) for f_name in files]
            
            for future in futures:
                future.result()
