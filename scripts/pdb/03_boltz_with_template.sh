#!/bin/bash
#SBATCH --job-name=pdb_boltz_with_templates
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=6-00:00:00
#SBATCH --output=logs/pdb_boltz_with_templates/slurm.%j.log

# env vars
export NXF_LOG_FILE=logs/pdb_boltz_with_templates/.nextflow.log
export NXF_CACHE_DIR=logs/pdb_boltz_with_templates/.nextflow

conda run -n nf-core --live-stream  nextflow run \
    ./workflows/pipelines/boltz.nf \
    --input "data/pdb/triad/staged/pdb_triad.cleaned.parquet" \
    -output-dir data/pdb \
    -profile gemini 