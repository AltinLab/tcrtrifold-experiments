#!/bin/bash
#SBATCH --job-name=extract_triad_conf_feat_test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=1-00:00:00
#SBATCH --output=tmp/nextflow/test/extract_triad_conf_feat.%j.log

if [ "${GITHUB_ACTIONS:-false}" = "true" ]; then
  PROFILE="gh_runner"
else
  PROFILE="standard"
fi

# env vars
export NXF_LOG_FILE=tmp/nextflow/test/extract_triad_conf_feat/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/test/extract_triad_conf_feat/cache

conda run -n nf-core --live-stream nextflow run \
    ./workflows/04_extract_feat.cresta.nf \
    --dset_name test \
    --inf_type af3