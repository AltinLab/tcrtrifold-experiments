#!/bin/bash
#SBATCH --job-name=pdb_tcrdock_boltz_template
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/pdb_tcrdock_boltz_template/slurm.%j.log

conda run -n tcrdock --live-stream python \
    ./workflows/subworkflows/local/extract_feat/resources/usr/bin/tcrdock_geom.py \
    --input_parquet data/pdb/triad/staged/pdb_triad.cleaned.parquet \
    --inference_type boltz \
    --topology_path data/pdb/triad/predictions_with_templates \
    --output_parquet_path data/pdb/triad/staged/pdb_triad.boltz_tcrdock_template.parquet