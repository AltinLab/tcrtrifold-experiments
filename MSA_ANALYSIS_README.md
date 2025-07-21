# MSA Analysis Pipeline

This pipeline evaluates the influence of paired and unpaired MSAs (Multiple Sequence Alignments) on confidence metrics and ability to distinguish binders and non-binders in TCR:p:MHC complexes.

## Overview

**Hypothesis**: Triads and p:MHCs that have higher Neff (number of effective sequences) and interface signal in their paired MSAs ("pairedMsa" in the AlphaFold MSA) have higher confidence metrics.

This is based on findings from:
- [AlphaFold2 complex prediction paper](https://www.nature.com/articles/s41467-022-28865-w)
- [AlphaFold3 paper](https://www.nature.com/articles/s41586-024-07487-w#Sec7)

## Pipeline Components

### 1. Core Analysis Module (`src/tcrtrifold/msa_analysis.py`)
- **Neff Calculation**: Implements the AlphaFold3 methodology for calculating effective sequence numbers
- **DCA Analysis**: Simplified Direct Coupling Analysis to measure interface signal strength
- **MSA Processing**: Functions to load and parse MSA data from JSON files
- **Interface Position Estimation**: Methods to identify interface residues

### 2. Processing Scripts
- **`workflows/bin/analyze_msa_metrics.py`**: Main script to process p:MHC and triad data
- **`workflows/bin/generate_msa_summary.py`**: Generates comprehensive summary statistics

### 3. Nextflow Pipeline (`workflows/07_msa_analysis.nf`)
Automated pipeline that:
- Processes both p:MHC complexes and full triads
- Extracts MSA metrics for all sequences
- Generates summary reports
- Supports sampling for testing

### 4. Analysis Notebook (`notebooks/msa_analysis_evaluation.ipynb`)
Jupyter notebook for:
- Data visualization and exploration
- Correlation analysis with confidence metrics
- Comparison of paired vs unpaired MSA effectiveness
- Statistical analysis and plots

## Usage

### Prerequisites

Activate the tcrtrifold-experiments conda environment:
```bash
conda env create -f envs/env_runner.yaml
conda activate tcrtrifold-experiments
```

### Quick Test

Run a quick test with sample data:
```bash
./scripts/test/test_msa_analysis.sh
```

### Full Pipeline (Nextflow)

Run the complete analysis:
```bash
# With sampling for testing
nextflow run workflows/07_msa_analysis.nf --sample_size 100

# Full dataset
nextflow run workflows/07_msa_analysis.nf
```

### Manual Analysis

Process specific files manually:
```bash
# Analyze p:MHC complexes
python workflows/bin/analyze_msa_metrics.py \
    --input data/test/pmhc/staged/test_pmhc.cleaned.parquet \
    --msa-dir data/test/msa \
    --output results/pmhc_msa_metrics.parquet \
    --data-type pmhc

# Analyze triads (TCR:p:MHC)
python workflows/bin/analyze_msa_metrics.py \
    --input data/test/triad/staged/test_triad.cleaned.parquet \
    --msa-dir data/test/msa \
    --output results/triad_msa_metrics.parquet \
    --data-type triad

# Generate summary
python workflows/bin/generate_msa_summary.py \
    --pmhc results/pmhc_msa_metrics.parquet \
    --triad results/triad_msa_metrics.parquet \
    --output results/msa_analysis_summary.json
```

## Data Structure

### Input Data
- **Parquet files**: Contain sequence data with columns:
  - `mhc_1_seq`, `mhc_2_seq`: MHC sequences
  - `tcr_1_seq`, `tcr_2_seq`: TCR sequences (triads only)
  - `peptide`: Peptide sequence
  - Confidence metrics (PAE, pLDDT, etc.)

- **MSA files**: JSON format in `data/test/msa/` containing:
  - `sequence`: Original sequence
  - `unpairedMsa`: Standard MSA in FASTA format
  - `pairedMsa`: Paired MSA in FASTA format
  - Located using SHA256 hash of sequences

### Output Data
- **MSA metrics parquet files**: Original data + MSA metrics
- **Summary JSON**: Statistical analysis and comparisons
- **Visualization plots**: From Jupyter notebook

## Metrics Calculated

### Neff (Number of Effective Sequences)
- Calculated using sequence identity threshold (80%)
- Applied to both paired and unpaired MSAs
- Based on AlphaFold3 methodology

### DCA Score (Interface Signal)
- Simplified Direct Coupling Analysis
- Measures covariation at interface positions
- Conservation-based scoring for interface residues

### Comparisons
- Paired vs unpaired MSA effectiveness
- Correlation with confidence metrics
- Statistical significance testing

## Key Outputs

1. **Neff Distributions**: Compare paired vs unpaired MSA depth
2. **Interface Signal**: DCA scores for binding interfaces
3. **Correlation Analysis**: MSA metrics vs confidence scores
4. **Statistical Comparisons**: Paired/unpaired ratios and significance

## Implementation Notes

- Uses polars for efficient data processing
- Handles missing MSA files gracefully
- Supports both MHC and TCR sequence types
- Interface positions estimated from sequence length (can be improved with structural data)
- Simplified DCA implementation (full DCA would require more sophisticated modeling)

## Future Improvements

1. **Full DCA Implementation**: Use proper statistical coupling analysis
2. **Structural Interface Definition**: Use PDB/AF3 structures to define true interface positions
3. **Binding Affinity Correlation**: Compare with experimental binding data
4. **Extended Confidence Metrics**: Analyze additional AF3 confidence scores
5. **Cross-validation**: Test on independent datasets

## References

- Evans, R. et al. Protein complex prediction with AlphaFold-Multimer. *Nature* 596, 590–596 (2022)
- Abramson, J. et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature* 630, 493–500 (2024)