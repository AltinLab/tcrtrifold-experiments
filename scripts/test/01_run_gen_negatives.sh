#!/bin/bash
#SBATCH --job-name=gen_negatives
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=8:00:00
#SBATCH --output=tmp/nextflow/test/gen_negatives.%j.log

if [ "${GITHUB_ACTIONS:-false}" = "true" ]; then
  PROFILE="gh_runner"
else
  PROFILE="standard"
fi

# env vars
export NXF_LOG_FILE=tmp/nextflow/test/gen_negatives/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/test/gen_negatives/cache

conda run -n nf-core --live-stream nextflow run \
    ./workflows/01_gen_negatives.nf \
        --dset_name test \
        --neg_depth 10 \
        --negs_from test \
        -profile "${PROFILE}"