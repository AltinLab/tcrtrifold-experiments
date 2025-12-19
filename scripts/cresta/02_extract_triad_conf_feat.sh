#!/bin/bash
#SBATCH --job-name=cresta_extract_triad_conf_feat
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/cresta_extract_triad_conf_feat/slurm.%j.log

conda run -n tcrtrifold-experiments --live-stream python \
    ./workflows/bin/extract_triad_conf_feat.py \
    --input_parquet data/cresta/triad/staged/cresta_triad.neg.parquet \
    --inference_type af3 \
    --inference_dir data/cresta/triad/inference \
    --output_path data/cresta/triad/cresta_triad.conf_af3.parquet