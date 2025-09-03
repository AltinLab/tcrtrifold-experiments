#!/bin/bash
#SBATCH --job-name=cresta
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/cresta/slurm.%j.log

# env vars
export NXF_LOG_FILE=logs/cresta/.nextflow.log
export NXF_CACHE_DIR=logs/cresta/.nextflow

conda run -n nf-core --live-stream  nextflow run \
    ./workflows/cresta.nf \
    --input data/cresta/raw/cresta_raw.csv \
    -output-dir data/cresta \
    -profile gemini \
    -resume