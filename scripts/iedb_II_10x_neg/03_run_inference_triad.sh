#!/bin/bash
#SBATCH --job-name=inference_triad
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=1-00:00:00
#SBATCH --output=tmp/nextflow/iedb_II_10x_neg/inference_triad.%j.log

# env vars
export NXF_LOG_FILE=tmp/nextflow/iedb_II_10x_neg/inference_triad/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/iedb_II_10x_neg/inference_triad/cache

conda run -n nf-core --live-stream nextflow run \
    ./workflows/03_inference_triad.nf \
        --dset_name iedb_II_10x_neg \
        --skip_msa 0 \
        --seeds 1,2,3,4,5