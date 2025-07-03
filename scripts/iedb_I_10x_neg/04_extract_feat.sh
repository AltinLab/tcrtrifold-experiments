#!/bin/bash
#SBATCH --job-name=extract_feat
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=1:00:00
#SBATCH --output=tmp/nextflow/iedb_I/extract_feat.%j.log

. ./scripts/setup.sh

# env vars
export NXF_LOG_FILE=tmp/nextflow/iedb_I/thresh_1.1x_neg/extract_feat/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/iedb_I/thresh_1.1x_neg/extract_feat/cache

# nextflow run \
#     ./workflows/04_extract_feat.iedb_I.thresh_1.1x_neg.nf \
#         -resume && \
nextflow run \
    ./workflows/04_recombine_feat.nf \
    --dset_name iedb_I \
    --input_pattern "*thresh_1*base*" \

# export NXF_LOG_FILE=tmp/nextflow/iedb_I/thresh_3.10x_neg/extract_feat/nextflow.log
# export NXF_CACHE_DIR=tmp/nextflow/iedb_I/thresh_3.10x_neg/extract_feat/cache

# nextflow run \
#     ./workflows/04_extract_feat.iedb_I.thresh_3.10x_neg.nf \
#         -resume && \
# nextflow run \
#     ./workflows/04_recombine_feat.nf \
#     --dset_name iedb_I \
#     --input_pattern "*thresh_3.10x_neg*base*" \