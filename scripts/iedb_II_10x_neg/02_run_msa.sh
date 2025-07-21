#!/bin/bash
#SBATCH --job-name=msa
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=1-00:00:00
#SBATCH --output=tmp/nextflow/iedb_II_10x_neg/msa.%j.log

# env vars
export NXF_LOG_FILE=tmp/nextflow/iedb_II_10x_neg/msa/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/iedb_II_10x_neg/msa/cache

conda run -n nf-core --live-stream nextflow run \
    ./workflows/02_msa.nf \
        --dset_name iedb_II_10x_neg