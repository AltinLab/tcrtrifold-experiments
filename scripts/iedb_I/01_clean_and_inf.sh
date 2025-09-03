#!/bin/bash
#SBATCH --job-name=iedb_I
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/iedb_I/slurm.%j.log

# env vars
export NXF_LOG_FILE=logs/iedb_I/.nextflow.log
export NXF_CACHE_DIR=logs/iedb_I/.nextflow

conda run -n nf-core --live-stream  nextflow run \
    ./workflows/iedb_I.nf \
    --input data/iedb_I/raw/immrep_IEDB.csv \
    --validation_exclusions data/pdb/triad/staged/pdb_validation_triad.annotated.parquet \
    -output-dir data/iedb_I \
    -profile gemini \
    -resume 