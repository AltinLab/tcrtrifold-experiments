#!/bin/bash
#SBATCH --job-name=iedb_II
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/iedb_II/slurm.%j.log

# env vars
export NXF_LOG_FILE=logs/iedb_II/.nextflow.log
export NXF_CACHE_DIR=logs/iedb_II/.nextflow

conda run -n nf-core --live-stream  nextflow run \
    ./workflows/iedb_II.nf \
    --input data/iedb_II/raw/immrep_IEDB.csv \
    --validation_exclusions data/cresta/triad/staged/cresta_triad.annotated.parquet \
    -output-dir data/iedb_II \
    -profile gemini \
    -resume