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
    ./workflows/subworkflows/local/cleaning/resources/usr/bin/clean_expansion.py \
    --clonotype_shortlist data/expansion/raw/75064_clonotypeShortlist.csv \
    --tcr_info data/expansion/raw/filtered_contig_annotations.csv \
    --hla_genotyping "data/expansion/raw/AHRI-UCSFtyped_HLA-genotyping_01-30-24 - Sheet1.csv" \
    --peptide_pool "data/expansion/raw/GigaPool_Sheet.xlsx - All Peptides.csv" \
    --output_triad_path data/expansion/triad/expansion_triad.parquet