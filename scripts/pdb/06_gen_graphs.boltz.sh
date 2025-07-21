#!/bin/bash
#SBATCH --job-name=gen_graphs_boltz_pdb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=1-00:00:00
#SBATCH --output=tmp/nextflow/pdb/gen_graphs_boltz.%j.log


# env vars
export NXF_LOG_FILE=tmp/nextflow/pdb/gen_graphs_boltz/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/pdb/gen_graphs_boltz/cache

conda run -n nf-core --live-stream nextflow run \
    ./workflows/06_gen_graphs.nf \
        --dset_name pdb \
        --inf_type boltz