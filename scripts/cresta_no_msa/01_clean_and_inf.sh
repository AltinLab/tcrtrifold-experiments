#!/bin/bash
#SBATCH --job-name=cresta_no_msa
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/cresta_no_msa/slurm.%j.log

# env vars
export NXF_LOG_FILE=logs/cresta_no_msa/.nextflow.log
export NXF_CACHE_DIR=logs/cresta_no_msa/.nextflow

conda run -n nf-core --live-stream  nextflow run \
    ./workflows/cresta_no_msa.nf \
    --input data/cresta/raw/cresta_raw.csv \
    -output-dir data/cresta_no_msa \
    -profile gemini \
    -resume