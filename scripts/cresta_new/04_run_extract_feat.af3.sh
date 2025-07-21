#!/bin/bash
#SBATCH --job-name=extract_triad_conf_feat_cresta_new
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=1-00:00:00
#SBATCH --output=tmp/nextflow/cresta_new/extract_triad_conf_feat.%j.log

# env vars
export NXF_LOG_FILE=tmp/nextflow/cresta_new/extract_triad_conf_feat/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/cresta_new/extract_triad_conf_feat/cache

conda run -n nf-core --live-stream nextflow run \
    ./workflows/04_extract_feat.cresta.nf \
    --dset_name cresta_new \
    --inf_type af3