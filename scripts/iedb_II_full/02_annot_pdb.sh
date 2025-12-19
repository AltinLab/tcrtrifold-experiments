#!/bin/bash
#SBATCH --job-name=iedb_II_full_annot_pdb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/iedb_II_full_annot_pdb/slurm.%j.log

# env vars
export NXF_LOG_FILE=logs/iedb_II_full_annot_pdb/.nextflow.log
export NXF_CACHE_DIR=logs/iedb_II_full_annot_pdb/.nextflow

conda run -n nf-core --live-stream  nextflow run \
    ./workflows/pipelines/annotate_pdb_matches.nf \
    --input data/iedb_II_full/triad/staged/iedb_II_triad.cleaned.parquet \
    -output-dir data/iedb_II_full \
    -profile gemini \
    -resume