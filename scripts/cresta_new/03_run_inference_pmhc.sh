#!/bin/bash
#SBATCH --job-name=inference_pmhc_cresta_new
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=5-00:00:00
#SBATCH --output=tmp/nextflow/cresta_new/inference_pmhc.%j.log

# env vars
export NXF_LOG_FILE=tmp/nextflow/cresta_new/inference_pmhc/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/cresta_new/inference_pmhc/cache

conda run -n nf-core --live-stream nextflow run \
    ./workflows/03_inference_pmhc.nf \
        --dset_name cresta_new \
        --skip_msa 0 \
        --seeds 1,2,3,4,5