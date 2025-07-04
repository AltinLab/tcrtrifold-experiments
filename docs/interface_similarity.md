# Interface Similarity Analysis

This module provides functionality to calculate interface similarity between protein complexes and their MSA alignments.

## Overview

The hypothesis is that triads (TCR:p:MHC complexes) and p:MHCs with higher "interface similarity" to their aligned paired sequences in the MSA (`pairedMsa`) may have different binding characteristics.

## Features

### Interface Similarity Metrics

For each protein chain in the complex, the following metrics are calculated:

- **mean_similarity**: Average similarity between the query interface sequence and all aligned sequences in the MSA
- **max_similarity**: Maximum similarity between the query interface sequence and any aligned sequence in the MSA  
- **num_sequences**: Number of aligned sequences used for comparison

### Interface Region Definitions

#### MHC Interface Regions
- Uses extracellular topological domains as defined in `tcrtrifold.mhc`
- Covers the regions that interact with TCRs and peptides
- Handles different HLA alleles (Class I and Class II)

#### TCR Interface Regions
- Uses CDR regions (CDR1, CDR2, CDR3) as defined by IMGT numbering
- These are the complementarity-determining regions that directly contact the MHC-peptide complex
- Extracted using `tcrtrifold.tcr` functionality

## Usage

### Python API

```python
from pathlib import Path
from tcrtrifold.interface_similarity import add_interface_similarity_features
import polars as pl

# Load data
df = pl.read_parquet("data/test/pmhc/staged/test_pmhc.cleaned.parquet")
msa_dir = Path("data/test/msa")

# Add interface similarity features
df_with_features = add_interface_similarity_features(
    df, msa_dir, complex_type="pmhc"
)

# For triads
triad_df = pl.read_parquet("data/test/triad/staged/test_triad.cleaned.parquet")
triad_with_features = add_interface_similarity_features(
    triad_df, msa_dir, complex_type="triad"
)
```

### Command Line

```bash
# Calculate interface similarity for a dataset
python workflows/bin/calculate_interface_similarity.py \
    --data_dir data \
    --dset_name test
```

### Nextflow Pipeline

```bash
# Run the interface similarity pipeline
nextflow run workflows/05_interface_similarity.nf \
    --data_dir data \
    --dset_name test \
    -profile gh_runner
```

## Output Columns

### For p:MHC complexes:
- `mhc_1_interface_mean_similarity`
- `mhc_1_interface_max_similarity`  
- `mhc_1_interface_num_sequences`
- `mhc_2_interface_mean_similarity`
- `mhc_2_interface_max_similarity`
- `mhc_2_interface_num_sequences`

### For triad complexes:
All of the above plus:
- `tcr_1_interface_mean_similarity`
- `tcr_1_interface_max_similarity`
- `tcr_1_interface_num_sequences`
- `tcr_2_interface_mean_similarity`
- `tcr_2_interface_max_similarity`
- `tcr_2_interface_num_sequences`

## Implementation Details

1. **Sequence Hashing**: Uses SHA256 hashing of sequences (via `tcrtrifold.utils.hash_sequence`) to locate corresponding MSA files
2. **MSA Loading**: Loads MSA data from JSON files in the format produced by the existing MSA pipeline
3. **Interface Extraction**: Extracts interface regions using existing topological domain definitions
4. **Similarity Calculation**: Uses `difflib.SequenceMatcher` to calculate sequence similarity ratios
5. **Batch Processing**: Handles large datasets efficiently by processing in batches

## Files

- `src/tcrtrifold/interface_similarity.py`: Main implementation
- `workflows/bin/calculate_interface_similarity.py`: Command-line script
- `workflows/05_interface_similarity.nf`: Nextflow pipeline
- `test_interface_similarity.py`: Test script

## Dependencies

- `polars`: For DataFrame operations
- `numpy`: For numerical computations
- `pathlib`: For file path handling
- `json`: For MSA file parsing
- `difflib`: For sequence similarity calculation
- `tcrtrifold.utils`: For sequence hashing
- `tcrtrifold.mhc`: For MHC interface extraction
- `tcrtrifold.tcr`: For TCR interface extraction