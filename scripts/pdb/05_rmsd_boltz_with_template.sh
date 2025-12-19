#!/bin/bash
#SBATCH --job-name=pdb_rmsd_boltz_template
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/pdb_rmsd_boltz_template/slurm.%j.log

conda run -n tcrtrifold-experiments --live-stream python \
    ./workflows/subworkflows/local/extract_feat/resources/usr/bin/rmsd.py \
    --input_parquet data/pdb/triad/staged/pdb_triad.cleaned.parquet \
    --cleaned_pdbs data/pdb/triad/cleaned_pdb \
    --inference_type boltz \
    --inference_dir data/pdb/triad/predictions_with_templates \
    --true_tcrdock_pq data/pdb/triad/staged/pdb_triad.true_tcrdock.parquet \
    --pred_tcrdock_pq data/pdb/triad/staged/pdb_triad.boltz_tcrdock_template.parquet \
    -o data/pdb/triad/staged/pdb_triad.boltz_template_rmsd.parquet