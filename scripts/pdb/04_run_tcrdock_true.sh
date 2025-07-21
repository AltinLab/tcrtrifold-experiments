#!/bin/bash
#SBATCH --job-name=tcrdock_true_pdb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=1-00:00:00
#SBATCH --output=tmp/nextflow/pdb/true_tcrdock.%j.log

# env vars
export NXF_LOG_FILE=tmp/nextflow/pdb/true_tcrdock/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/pdb/true_tcrdock/cache

conda run -n nf-core --live-stream nextflow run \
    ./workflows/04_tcrdock.nf \
    --dset_name pdb \
    --from_true_struct