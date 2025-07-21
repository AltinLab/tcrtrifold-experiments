#!/bin/bash
#SBATCH --job-name=rmsd_boltz_pdb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=1-00:00:00
#SBATCH --output=tmp/nextflow/pdb/rmsd_boltz.%j.log

# env vars
export NXF_LOG_FILE=tmp/nextflow/pdb/rmsd_boltz/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/pdb/rmsd_boltz/cache

conda run -n nf-core --live-stream nextflow run \
    ./workflows/05_rmsd.nf \
    --inf_type boltz