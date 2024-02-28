import os
import subprocess
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

def process_row(row, src_base_path, src):
    original_name = row['original_name']
    genome_file_name = row['genome_file_name']
    gff_file_name = row['gff_file_name']

    output_file = f'data/{src}/cds_from_gff/'
    if not os.path.exists(output_file):
        os.makedirs(output_file)

    if os.path.exists(f'{src_base_path}/cds_from_gff/{original_name}_cds_from_gff.fna'):
        return

    gffread_command = ['src/utils/gffread/gffread/gffread', '-x', f'{src_base_path}/cds_from_gff/{original_name}_cds_from_gff.fna', '-g', f'{src_base_path}/genomes/{genome_file_name}', f'{src_base_path}/gff/{gff_file_name}', '-W', '--attrs', 'CDS']
    process = subprocess.Popen(gffread_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, stderr = process.communicate()

    if stderr:
        err = stderr.decode()

        if "Error" in err:
            print(f'Encountered an error: {err}')

            if os.path.exists(f'{src_base_path}/genomes/{genome_file_name}.fai'):
                os.remove(f'{src_base_path}/genomes/{genome_file_name}.fai')
            return

        if "Warning" in err:
            print(original_name, f'{src_base_path}/gff/{gff_file_name}')

    if os.path.exists(f'{src_base_path}/genomes/{genome_file_name}.fai'):
        os.remove(f'{src_base_path}/genomes/{genome_file_name}.fai')


def process_rows_chunk(df_chunk, src_base_path, src):
    for _, row in df_chunk.iterrows():
        process_row(row, src_base_path, src)


def create_cds_from_gff():
    sources = ['NCBI', 'FungiDB', 'EnsemblFungi', 'MycoCosm']
    base_paths = ['data/NCBI', 'data/FungiDB', 'data/EnsemblFungi', 'data/MycoCosm']
    manifests = ['ncbi_input_species.csv', 'fungidb_input_species.csv', 'ensemblfungi_input_species.csv', 'mycocosm_input_species.csv']

    num_cores = max(1, os.cpu_count() - 1) # Use all cores except one
    for m_idx, m in enumerate(manifests):
        count = 0
        csv_path = os.path.join(base_paths[m_idx], m)

        df = pd.read_csv(csv_path).drop_duplicates(subset='cds_file_name', ignore_index=True)
        
        chunk_size = max(1, len(df) // num_cores) # Calculate chunk size based on number of cores
        chunks = [df[i:i + chunk_size] for i in range(0, df.shape[0], chunk_size)]

        with ThreadPoolExecutor(max_workers=num_cores) as executor:
            futures = []
            for chunk in chunks:
                futures.append(executor.submit(process_rows_chunk, chunk, base_paths[m_idx], sources[m_idx]))

            for _ in as_completed(futures):
                count += chunk_size
                print(f'Done with {min(count, df.shape[0])} {sources[m_idx]} species...')

        files = os.listdir(f'data/{sources[m_idx]}/cds_from_gff')
        df['cds_file_name'] = df['cds_file_name'].str.replace('_cds', '_cds_from_gff')

        df = df[df['cds_file_name'].isin(files)]

        df.to_csv(f'data/{sources[m_idx]}/{sources[m_idx].lower()}_gff_input_species.csv', index=False)
    print()