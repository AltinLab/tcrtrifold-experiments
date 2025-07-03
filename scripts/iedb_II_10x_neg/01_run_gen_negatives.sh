#!/bin/bash
#SBATCH --job-name=gen_negatives
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=1-00:00:00
#SBATCH --output=tmp/nextflow/iedb_II_10x_neg/gen_negatives.%j.log

# env vars
export NXF_LOG_FILE=tmp/nextflow/iedb_II_10x_neg/gen_negatives/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/iedb_II_10x_neg/gen_negatives/cache

conda run -n nf-core --live-stream nextflow run \
    ./workflows/01_gen_negatives.nf \
        --dset_name iedb_II_10x_neg \
        --neg_depth 10 \
        --negs_from iedb_II_1x_neg,iedb_I_1x_neg