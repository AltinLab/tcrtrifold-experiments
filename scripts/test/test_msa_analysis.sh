#!/bin/bash
# Test script for MSA analysis pipeline
# 
# This script demonstrates how to run the MSA analysis pipeline
# on the test dataset.

set -e

echo "MSA Analysis Pipeline Test"
echo "========================="

# Check if we're in the correct directory
if [ ! -f "workflows/07_msa_analysis.nf" ]; then
    echo "Error: Please run this script from the tcrtrifold-experiments root directory"
    exit 1
fi

# Check if conda environment is available
if ! conda env list | grep -q "tcrtrifold-experiments"; then
    echo "Warning: tcrtrifold-experiments conda environment not found"
    echo "Please create it with: conda env create -f envs/env_runner.yaml"
fi

echo "Testing with sample data (first 10 rows)..."

# Create results directory
mkdir -p results/msa_analysis_test

# Test p:MHC analysis with small sample
echo "1. Testing p:MHC analysis..."
python workflows/bin/analyze_msa_metrics.py \
    --input data/test/pmhc/staged/test_pmhc.cleaned.parquet \
    --msa-dir data/test/msa \
    --output results/msa_analysis_test/pmhc_test.parquet \
    --data-type pmhc \
    --sample-size 10

# Test triad analysis with small sample
echo "2. Testing triad analysis..."
python workflows/bin/analyze_msa_metrics.py \
    --input data/test/triad/staged/test_triad.cleaned.parquet \
    --msa-dir data/test/msa \
    --output results/msa_analysis_test/triad_test.parquet \
    --data-type triad \
    --sample-size 10

# Generate summary
echo "3. Generating summary..."
python workflows/bin/generate_msa_summary.py \
    --pmhc results/msa_analysis_test/pmhc_test.parquet \
    --triad results/msa_analysis_test/triad_test.parquet \
    --output results/msa_analysis_test/summary.json

echo "Test completed successfully!"
echo "Results in: results/msa_analysis_test/"
echo ""
echo "To run the full pipeline with Nextflow:"
echo "  nextflow run workflows/07_msa_analysis.nf --sample_size 100"
echo ""
echo "To run without sampling (full dataset):"
echo "  nextflow run workflows/07_msa_analysis.nf"