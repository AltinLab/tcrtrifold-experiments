#!/bin/bash
#SBATCH --job-name=msa
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=1-00:00:00
#SBATCH --output=tmp/nextflow/test/msa.%j.log

if [ "${GITHUB_ACTIONS:-false}" = "true" ]; then
  PROFILE="gh_runner"
else
  PROFILE="standard"
fi

# env vars
export NXF_LOG_FILE=tmp/nextflow/test/msa/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/test/msa/cache

conda run -n nf-core --live-stream nextflow run \
    ./workflows/02_msa.nf \
        --dset_name test