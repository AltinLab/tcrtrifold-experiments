#!/bin/bash
#SBATCH --job-name=gen_negatives
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=1-00:00:00
#SBATCH --output=tmp/nextflow/pdb/gen_negatives.%j.log

# env vars
export NXF_LOG_FILE=tmp/nextflow/pdb/gen_negatives/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/pdb/gen_negatives/cache

conda run -n nf-core --live-stream nextflow run \
    ./workflows/01_gen_negatives.nf \
        --dset_name pdb \
        --neg_depth 10 \
        --negs_from pdb