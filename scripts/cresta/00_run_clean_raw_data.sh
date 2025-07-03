#!/bin/bash
#SBATCH --job-name=clean_raw_data
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=1:00:00
#SBATCH --output=tmp/nextflow/cresta/clean_raw_data.%j.log

. ./scripts/setup.sh

# env vars
export NXF_LOG_FILE=tmp/nextflow/cresta/clean_raw_data/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/cresta/clean_raw_data/cache

nextflow run \
    ./workflows/00_clean_raw_data.cresta.nf