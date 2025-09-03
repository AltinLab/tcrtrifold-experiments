#!/bin/bash
#SBATCH --job-name=iedb_I_extract_triad_conf_feat
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/iedb_I_extract_triad_conf_feat/slurm.%j.log

conda run -n tcrtrifold-experiments --live-stream python \
    ./workflows/bin/extract_triad_conf_feat.py \
    --input_parquet data/iedb_I/triad/staged/iedb_I_triad.neg.parquet \
    --inference_type af3 \
    --inference_dir data/iedb_I/triad/inference \
    --output_path data/iedb_I/triad/iedb_I_triad.conf_af3.parquet