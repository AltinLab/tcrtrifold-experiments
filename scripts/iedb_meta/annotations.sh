#!/bin/bash
#SBATCH --job-name=iedb_meta_annot
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/iedb_meta_annot/slurm.%j.log

# env vars
export NXF_LOG_FILE=logs/iedb_meta_annot/.nextflow.log
export NXF_CACHE_DIR=logs/iedb_meta_annot/.nextflow

conda run -n nf-core --live-stream  nextflow run \
    ./workflows/pipelines/annotate_pdb_matches.nf \
    --input data/iedb_meta/addtl_xray_corrected.parquet \
    -output-dir data/iedb_meta \
    -profile gemini