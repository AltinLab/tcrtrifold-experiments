#!/usr/bin/env python3
"""
Script to analyze MSA metrics for p:MHC and triad data.

This script processes parquet files containing sequence data and extracts
MSA-based metrics including Neff and interface signal for both paired
and unpaired MSAs.
"""

import argparse
import sys
from pathlib import Path
import polars as pl

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parents[2] / "src"))

from tcrtrifold.msa_analysis import analyze_msa_metrics, estimate_interface_positions


def process_pmhc_data(
    input_file: Path,
    msa_dir: Path,
    output_file: Path,
    sample_size: int = None
) -> None:
    """
    Process p:MHC data and extract MSA metrics.
    
    Args:
        input_file (Path): Input parquet file with p:MHC data
        msa_dir (Path): Directory containing MSA files
        output_file (Path): Output file for results
        sample_size (int, optional): Sample size for testing
    """
    print(f"Processing p:MHC data from {input_file}")
    
    # Load data
    df = pl.read_parquet(input_file)
    
    if sample_size:
        df = df.head(sample_size)
        print(f"Sampling {sample_size} rows for testing")
    
    print(f"Loaded {len(df)} p:MHC complexes")
    
    results = []
    
    for i, row in enumerate(df.iter_rows(named=True)):
        if i % 100 == 0:
            print(f"Processing row {i+1}/{len(df)}")
        
        # Extract sequences
        peptide = row.get('peptide', '')
        mhc_1_seq = row.get('mhc_1_seq', '')
        mhc_2_seq = row.get('mhc_2_seq', '')
        
        # Analyze MHC sequences
        mhc_1_metrics = {}
        mhc_2_metrics = {}
        
        if mhc_1_seq:
            interface_pos = estimate_interface_positions(len(mhc_1_seq))
            mhc_1_metrics = analyze_msa_metrics(
                mhc_1_seq, msa_dir, seq_type='mhc', interface_positions=interface_pos
            )
            # Add prefix to distinguish metrics
            mhc_1_metrics = {f"mhc_1_{k}": v for k, v in mhc_1_metrics.items()}
        
        if mhc_2_seq:
            interface_pos = estimate_interface_positions(len(mhc_2_seq))
            mhc_2_metrics = analyze_msa_metrics(
                mhc_2_seq, msa_dir, seq_type='mhc', interface_positions=interface_pos
            )
            # Add prefix to distinguish metrics
            mhc_2_metrics = {f"mhc_2_{k}": v for k, v in mhc_2_metrics.items()}
        
        # Combine all data
        result_row = {
            **row,
            **mhc_1_metrics,
            **mhc_2_metrics,
            'peptide_length': len(peptide),
            'mhc_1_length': len(mhc_1_seq),
            'mhc_2_length': len(mhc_2_seq),
        }
        
        results.append(result_row)
    
    # Create output DataFrame
    result_df = pl.DataFrame(results)
    
    # Save results
    result_df.write_parquet(output_file)
    print(f"Results saved to {output_file}")
    
    # Print summary statistics
    print("\nSummary Statistics:")
    print(f"Total complexes processed: {len(result_df)}")
    
    for col in ['mhc_1_unpaired_neff', 'mhc_1_paired_neff', 'mhc_2_unpaired_neff', 'mhc_2_paired_neff']:
        if col in result_df.columns:
            stats = result_df.select(pl.col(col)).describe()
            print(f"{col}: mean={stats['mean'][0]:.2f}, std={stats['std'][0]:.2f}")


def process_triad_data(
    input_file: Path,
    msa_dir: Path,
    output_file: Path,
    sample_size: int = None
) -> None:
    """
    Process triad data and extract MSA metrics.
    
    Args:
        input_file (Path): Input parquet file with triad data
        msa_dir (Path): Directory containing MSA files
        output_file (Path): Output file for results
        sample_size (int, optional): Sample size for testing
    """
    print(f"Processing triad data from {input_file}")
    
    # Load data
    df = pl.read_parquet(input_file)
    
    if sample_size:
        df = df.head(sample_size)
        print(f"Sampling {sample_size} rows for testing")
    
    print(f"Loaded {len(df)} triads")
    
    results = []
    
    for i, row in enumerate(df.iter_rows(named=True)):
        if i % 100 == 0:
            print(f"Processing row {i+1}/{len(df)}")
        
        # Extract sequences
        peptide = row.get('peptide', '')
        mhc_1_seq = row.get('mhc_1_seq', '')
        mhc_2_seq = row.get('mhc_2_seq', '')
        tcr_1_seq = row.get('tcr_1_seq', '')
        tcr_2_seq = row.get('tcr_2_seq', '')
        
        # Analyze MHC sequences
        mhc_1_metrics = {}
        mhc_2_metrics = {}
        
        if mhc_1_seq:
            interface_pos = estimate_interface_positions(len(mhc_1_seq))
            mhc_1_metrics = analyze_msa_metrics(
                mhc_1_seq, msa_dir, seq_type='mhc', interface_positions=interface_pos
            )
            mhc_1_metrics = {f"mhc_1_{k}": v for k, v in mhc_1_metrics.items()}
        
        if mhc_2_seq:
            interface_pos = estimate_interface_positions(len(mhc_2_seq))
            mhc_2_metrics = analyze_msa_metrics(
                mhc_2_seq, msa_dir, seq_type='mhc', interface_positions=interface_pos
            )
            mhc_2_metrics = {f"mhc_2_{k}": v for k, v in mhc_2_metrics.items()}
        
        # Analyze TCR sequences
        tcr_1_metrics = {}
        tcr_2_metrics = {}
        
        if tcr_1_seq:
            interface_pos = estimate_interface_positions(len(tcr_1_seq))
            tcr_1_metrics = analyze_msa_metrics(
                tcr_1_seq, msa_dir, seq_type='tcr', interface_positions=interface_pos
            )
            tcr_1_metrics = {f"tcr_1_{k}": v for k, v in tcr_1_metrics.items()}
        
        if tcr_2_seq:
            interface_pos = estimate_interface_positions(len(tcr_2_seq))
            tcr_2_metrics = analyze_msa_metrics(
                tcr_2_seq, msa_dir, seq_type='tcr', interface_positions=interface_pos
            )
            tcr_2_metrics = {f"tcr_2_{k}": v for k, v in tcr_2_metrics.items()}
        
        # Combine all data
        result_row = {
            **row,
            **mhc_1_metrics,
            **mhc_2_metrics,
            **tcr_1_metrics,
            **tcr_2_metrics,
            'peptide_length': len(peptide),
            'mhc_1_length': len(mhc_1_seq),
            'mhc_2_length': len(mhc_2_seq),
            'tcr_1_length': len(tcr_1_seq),
            'tcr_2_length': len(tcr_2_seq),
        }
        
        results.append(result_row)
    
    # Create output DataFrame
    result_df = pl.DataFrame(results)
    
    # Save results
    result_df.write_parquet(output_file)
    print(f"Results saved to {output_file}")
    
    # Print summary statistics
    print("\nSummary Statistics:")
    print(f"Total triads processed: {len(result_df)}")
    
    for prefix in ['mhc_1', 'mhc_2', 'tcr_1', 'tcr_2']:
        for suffix in ['unpaired_neff', 'paired_neff']:
            col = f"{prefix}_{suffix}"
            if col in result_df.columns:
                stats = result_df.select(pl.col(col)).describe()
                print(f"{col}: mean={stats['mean'][0]:.2f}, std={stats['std'][0]:.2f}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze MSA metrics for p:MHC and triad data"
    )
    parser.add_argument(
        "--input", "-i",
        type=Path,
        required=True,
        help="Input parquet file"
    )
    parser.add_argument(
        "--msa-dir", "-m",
        type=Path,
        required=True,
        help="Directory containing MSA files"
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        required=True,
        help="Output parquet file"
    )
    parser.add_argument(
        "--data-type", "-t",
        choices=['pmhc', 'triad'],
        required=True,
        help="Type of data to process"
    )
    parser.add_argument(
        "--sample-size", "-s",
        type=int,
        default=None,
        help="Sample size for testing (use all data if not specified)"
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    if not args.input.exists():
        print(f"Error: Input file {args.input} does not exist")
        sys.exit(1)
    
    if not args.msa_dir.exists():
        print(f"Error: MSA directory {args.msa_dir} does not exist")
        sys.exit(1)
    
    # Create output directory if needed
    args.output.parent.mkdir(parents=True, exist_ok=True)
    
    # Process data based on type
    if args.data_type == 'pmhc':
        process_pmhc_data(args.input, args.msa_dir, args.output, args.sample_size)
    elif args.data_type == 'triad':
        process_triad_data(args.input, args.msa_dir, args.output, args.sample_size)


if __name__ == "__main__":
    main()