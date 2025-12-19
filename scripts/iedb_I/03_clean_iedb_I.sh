#!/bin/bash
#SBATCH --job-name=iedb_I_receptor_ref
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/iedb_I_receptor_ref/slurm.%j.log

conda run -n tcrtrifold-experiments --live-stream python \
    workflows/subworkflows/local/cleaning/resources/usr/bin/clean_iedb_I.py \
    -c data/iedb_I/raw/immrep_IEDB.csv \
    -op /dev/null \
    -ot data/iedb_I/triad/iedb_I_triad.receptor_reference.parquet