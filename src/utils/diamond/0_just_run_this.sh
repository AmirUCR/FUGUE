#!/bin/bash
conda run -n fugue python 1_create_sh_script.py
conda run -n fugue python 2_make_proteome.py
conda run -n diamond2 ./3_makedb_for_ref_species.sh
conda run -n diamond2 ./4_makedb_for_all_species.sh
conda run -n diamond2 ./5_run_bidirectional_diamond.sh
conda run -n fugue python 6_aggregate_orthogroups.py