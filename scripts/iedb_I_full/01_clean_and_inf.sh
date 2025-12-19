#!/bin/bash
#SBATCH --job-name=iedb_I
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/iedb_I_full/slurm.%j.log

# env vars
export NXF_LOG_FILE=logs/iedb_I_full/.nextflow.log
export NXF_CACHE_DIR=logs/iedb_I_full/.nextflow

conda run -n nf-core --live-stream  nextflow run \
    ./workflows/iedb_I_full.nf \
    --input data/iedb_I/raw/immrep_IEDB.csv \
    -output-dir data/iedb_I_full \
    -profile gemini \
    -resume 