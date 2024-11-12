#!/bin/bash
# Generated script by 1_create_sh_script.py

PROT_PATH="/home/amohs002/projects/research/ALLEGRO_Fungi_Downloader/data/fourdbs_concat/proteomes/saccharomyces_cerevisiae.faa"
FILENAME=saccharomyces_cerevisiae

diamond makedb --in $PROT_PATH -d inputs/reference/$FILENAME
