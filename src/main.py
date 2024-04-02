import os
import sys
import signal
import threading
import subprocess

from utils.merger import merge_dbs, merge_gffs
from utils.delimited_ortho import delimited_ortho
from utils.cds_from_gff import create_cds_from_gff
from utils.cds_from_gff_delimiter import cds_from_gff_delimiter
from downloaders.NCBI.ncbi_downloader import NCBI_Downloader
from downloaders.FungiDB.fungidb_downloader import FungiDB_Downloader
from downloaders.EnsemblFungi.ensembl_download import EnsemblFungi_Downloader
from downloaders.MycoCosm.mycocosm_download import MycoCosm_Downloader


def signal_handler(sig, frame):
    print(f'\nExiting FUGUE. Ciao.')
    sys.exit(0)


def run_diamond():
    script_dir = 'src/utils/diamond/'
    script_path = '0_just_run_this.sh'

    # Open a subprocess to run the shell script within its own directory
    process = subprocess.Popen(['bash', script_path], cwd=script_dir)
    process.wait()


def run_ortho_finder():
    script_dir = 'src/utils/ortholog_finder/'
    script_path = 'find_orthogroup.py'

    # Open a subprocess to run the shell script within its own directory
    process = subprocess.Popen(['python', script_path], cwd=script_dir)
    process.wait()


def initialize_downloader(downloader_class):
    downloader = downloader_class()
    downloader.download()


def main(choice_arg=''):
    signal.signal(signal.SIGINT, signal_handler)

    if not os.path.exists('data'):
        os.mkdir('data')

    choice = '0'

    downloaders = {
        '1': NCBI_Downloader,
        '2': FungiDB_Downloader,
        '3': EnsemblFungi_Downloader,
        '4': MycoCosm_Downloader,
    }

    if choice_arg:
        print(f'Received choice {choice_arg} via argv.')
        choices = choice_arg
    else:
        print('Welcome to FUGUE. Select an option:\n')
        choices = input('1. NCBI Datasets\n2. FungiDB\n3. EnsemblFungi\n4. MycoCosm\n5. Download all\n6. Create CDS from GFF\n7. Create Delimited CDS from GFF\n8. Merge Downloads\n9. Merge GFF\n10. Run DIAMOND\n11. Run Orthofinder\n12. Delimited Ortho\n13. All of the Above\n14. Quit\nYour choice(s) (separate by &, for example: 6 & 7 & 8 & 9): ')

    choices = choices.split(' & ')

    for choice in choices:
        if choice in downloaders:
            downloader = downloaders[choice]()
            downloader.download()

        if choice == '5':
            threads = []
            for _, downloader_class in downloaders.items():
                thread = threading.Thread(target=initialize_downloader, args=(downloader_class,))
                threads.append(thread)
                thread.start()

            for thread in threads:
                thread.join()

        elif choice == '6':
            create_cds_from_gff()

        elif choice == '7':
            cds_from_gff_delimiter()

        elif choice == '8':
            merge_dbs()

        elif choice == '9':
            merge_gffs()
        
        elif choice == '10':
            run_diamond()

        elif choice == '11':
            run_ortho_finder()

        elif choice == '12':
            delimited_ortho()

        elif choice == '13':
            threads = []
            for _, downloader_class in downloaders.items():
                thread = threading.Thread(target=initialize_downloader, args=(downloader_class,))
                threads.append(thread)
                thread.start()

            for thread in threads:
                thread.join()
            
            print('Creating CDS from GFF...')
            create_cds_from_gff()

            print('Creating delimited CDS from GFF...')
            cds_from_gff_delimiter()

            print('Merging DBs...')
            merge_dbs()

            print('Merging GFFs')
            merge_gffs()

            print('DIAMOND')
            run_diamond()

            print('Orthofinder')
            run_ortho_finder()

            print('Creating delimited ortho')
            delimited_ortho()

    print('Goodbye.')
    return 0
    

if __name__ == '__main__':
    if len(sys.argv) > 1:
        sys.exit(main(sys.argv[1]))
    sys.exit(main())