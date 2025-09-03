#!/bin/bash
#SBATCH --job-name=pdb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/pdb/slurm.%j.log

# env vars
export NXF_LOG_FILE=logs/pdb/.nextflow.log
export NXF_CACHE_DIR=logs/pdb/.nextflow

conda run -n nf-core --live-stream  nextflow run \
    ./workflows/pdb.nf \
    --input_replication "data/pdb/raw/table_S1_structure_benchmark_complexes.csv" \
    --input_stcr "data/pdb/raw/db_summary.dat" \
    -output-dir data/pdb \
    -profile gemini \
    -resume