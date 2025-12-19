#!/bin/bash
#SBATCH --job-name=iedb_II_corrected
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/iedb_II_corrected/slurm.%j.log

# env vars
export NXF_LOG_FILE=logs/iedb_II_corrected/.nextflow.log
export NXF_CACHE_DIR=logs/iedb_II_corrected/.nextflow

conda run -n nf-core --live-stream  nextflow run \
    ./workflows/iedb_II_corrected.nf \
    --input_triad data/iedb_II_corrected/triad/iedb_II_corrected_triad.cleaned.parquet \
    --input_antigen data/iedb_II_corrected/triad/iedb_II_corrected_triad.cleaned.parquet \
    --validation_exclusions data/cresta/triad/staged/cresta_triad.annotated.parquet \
    -output-dir data/iedb_II_corrected \
    -profile gemini