#!/bin/bash
#SBATCH --job-name=clean_raw_data
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=1:00:00
#SBATCH --output=tmp/nextflow/cresta_new/clean_raw_data.%j.log

# env vars
export NXF_LOG_FILE=tmp/nextflow/cresta_new/clean_raw_data/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/cresta_new/clean_raw_data/cache

conda run -n nf-core --live-stream  nextflow run \
    ./workflows/00_clean_raw_data.cresta.nf \
    --dset_name cresta_new