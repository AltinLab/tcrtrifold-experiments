#!/bin/bash
#SBATCH --job-name=interface_similarity
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=1-00:00:00
#SBATCH --output=tmp/nextflow/test/interface_similarity.%j.log

# Since nextflow plugins have network issues in this environment,
# run the interface similarity calculation directly with Python

mkdir -p tmp/nextflow/test/interface_similarity
export PYTHONPATH="${PWD}/src:${PYTHONPATH}"

echo "Running interface similarity analysis directly with Python..."
conda run -n tcrtrifold-experiments python workflows/bin/calculate_interface_similarity.py \
    --data_dir data \
    --dset_name test

echo "Interface similarity analysis completed successfully!"