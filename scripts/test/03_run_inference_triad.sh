#!/bin/bash
#SBATCH --job-name=inference_triad
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=1-00:00:00
#SBATCH --output=tmp/nextflow/test/inference_triad.%j.log

if [ "${GITHUB_ACTIONS:-false}" = "true" ]; then
  PROFILE="gh_runner"
else
  PROFILE="standard"
fi

# env vars
export NXF_LOG_FILE=tmp/nextflow/test/inference_triad/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/test/inference_triad/cache

conda run -n nf-core --live-stream nextflow run \
    ./workflows/03_inference_triad.nf \
        --dset_name test \
        --skip_msa 0