#!/bin/bash
#SBATCH --job-name=iedb_II
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/iedb_II_full/slurm.%j.log

# env vars
export NXF_LOG_FILE=logs/iedb_II_full/.nextflow.log
export NXF_CACHE_DIR=logs/iedb_II_full/.nextflow

conda run -n nf-core --live-stream  nextflow run \
    ./workflows/iedb_II_full.nf \
    --input data/iedb_II/raw/immrep_IEDB.csv \
    -output-dir data/iedb_II_full \
    -profile gemini \
    -resume