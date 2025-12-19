#!/bin/bash
#SBATCH --job-name=pdb_extract_triad_conf_feat
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/pdb_extract_triad_conf_feat/slurm.%j.log

conda run -n tcrtrifold-experiments --live-stream python \
    ./workflows/bin/extract_triad_conf_feat.py \
    --input_parquet data/pdb/triad/staged/pdb_triad.cleaned.parquet \
    --inference_type af3 \
    --inference_dir data/pdb/triad/inference \
    --summary_only \
    --output_path data/pdb/triad/pdb_triad.conf_af3.parquet