#!/usr/bin/env python3
"""
Script to generate summary statistics for MSA analysis results.
"""

import argparse
import json
import sys
from pathlib import Path
import polars as pl
import numpy as np

def calculate_summary_stats(df: pl.DataFrame, prefix: str) -> dict:
    """Calculate summary statistics for MSA metrics."""
    stats = {}
    
    # Metrics to analyze
    metrics = [
        'unpaired_neff', 'paired_neff', 
        'unpaired_dca_score', 'paired_dca_score'
    ]
    
    # Sequence types to analyze
    seq_types = []
    if prefix == 'pmhc':
        seq_types = ['mhc_1', 'mhc_2']
    elif prefix == 'triad':
        seq_types = ['mhc_1', 'mhc_2', 'tcr_1', 'tcr_2']
    
    for seq_type in seq_types:
        for metric in metrics:
            col_name = f"{seq_type}_{metric}"
            if col_name in df.columns:
                values = df.select(pl.col(col_name)).to_numpy().flatten()
                values = values[~np.isnan(values)]  # Remove NaN values
                
                if len(values) > 0:
                    stats[col_name] = {
                        'count': len(values),
                        'mean': float(np.mean(values)),
                        'std': float(np.std(values)),
                        'min': float(np.min(values)),
                        'max': float(np.max(values)),
                        'median': float(np.median(values)),
                        'q25': float(np.percentile(values, 25)),
                        'q75': float(np.percentile(values, 75))
                    }
    
    return stats

def compare_paired_vs_unpaired(df: pl.DataFrame, prefix: str) -> dict:
    """Compare paired vs unpaired MSA metrics."""
    comparisons = {}
    
    seq_types = []
    if prefix == 'pmhc':
        seq_types = ['mhc_1', 'mhc_2']
    elif prefix == 'triad':
        seq_types = ['mhc_1', 'mhc_2', 'tcr_1', 'tcr_2']
    
    for seq_type in seq_types:
        paired_neff_col = f"{seq_type}_paired_neff"
        unpaired_neff_col = f"{seq_type}_unpaired_neff"
        paired_dca_col = f"{seq_type}_paired_dca_score"
        unpaired_dca_col = f"{seq_type}_unpaired_dca_score"
        
        if all(col in df.columns for col in [paired_neff_col, unpaired_neff_col]):
            # Compare Neff values
            paired_neff = df.select(pl.col(paired_neff_col)).to_numpy().flatten()
            unpaired_neff = df.select(pl.col(unpaired_neff_col)).to_numpy().flatten()
            
            # Remove NaN values
            valid_mask = ~(np.isnan(paired_neff) | np.isnan(unpaired_neff))
            paired_neff = paired_neff[valid_mask]
            unpaired_neff = unpaired_neff[valid_mask]
            
            if len(paired_neff) > 0:
                neff_ratio = paired_neff / (unpaired_neff + 1e-10)  # Avoid division by zero
                comparisons[f"{seq_type}_neff_comparison"] = {
                    'paired_higher_count': int(np.sum(paired_neff > unpaired_neff)),
                    'unpaired_higher_count': int(np.sum(unpaired_neff > paired_neff)),
                    'equal_count': int(np.sum(paired_neff == unpaired_neff)),
                    'mean_ratio_paired_to_unpaired': float(np.mean(neff_ratio)),
                    'median_ratio_paired_to_unpaired': float(np.median(neff_ratio))
                }
        
        if all(col in df.columns for col in [paired_dca_col, unpaired_dca_col]):
            # Compare DCA scores
            paired_dca = df.select(pl.col(paired_dca_col)).to_numpy().flatten()
            unpaired_dca = df.select(pl.col(unpaired_dca_col)).to_numpy().flatten()
            
            # Remove NaN values
            valid_mask = ~(np.isnan(paired_dca) | np.isnan(unpaired_dca))
            paired_dca = paired_dca[valid_mask]
            unpaired_dca = unpaired_dca[valid_mask]
            
            if len(paired_dca) > 0:
                dca_ratio = paired_dca / (unpaired_dca + 1e-10)
                comparisons[f"{seq_type}_dca_comparison"] = {
                    'paired_higher_count': int(np.sum(paired_dca > unpaired_dca)),
                    'unpaired_higher_count': int(np.sum(unpaired_dca > paired_dca)),
                    'equal_count': int(np.sum(paired_dca == unpaired_dca)),
                    'mean_ratio_paired_to_unpaired': float(np.mean(dca_ratio)),
                    'median_ratio_paired_to_unpaired': float(np.median(dca_ratio))
                }
    
    return comparisons

def analyze_confidence_correlation(df: pl.DataFrame) -> dict:
    """Analyze correlation between MSA metrics and confidence metrics."""
    correlations = {}
    
    # Look for confidence metrics in the data
    confidence_cols = []
    for col in df.columns:
        if any(keyword in col.lower() for keyword in ['pae', 'confidence', 'lddt', 'plddt']):
            confidence_cols.append(col)
    
    if not confidence_cols:
        return {'message': 'No confidence metrics found in data'}
    
    # MSA metric columns
    msa_cols = []
    for col in df.columns:
        if any(keyword in col for keyword in ['neff', 'dca_score']):
            msa_cols.append(col)
    
    for conf_col in confidence_cols:
        for msa_col in msa_cols:
            try:
                # Calculate correlation
                conf_vals = df.select(pl.col(conf_col)).to_numpy().flatten()
                msa_vals = df.select(pl.col(msa_col)).to_numpy().flatten()
                
                # Remove NaN values
                valid_mask = ~(np.isnan(conf_vals) | np.isnan(msa_vals))
                conf_vals = conf_vals[valid_mask]
                msa_vals = msa_vals[valid_mask]
                
                if len(conf_vals) > 10:  # Need sufficient data points
                    correlation = np.corrcoef(conf_vals, msa_vals)[0, 1]
                    correlations[f"{conf_col}_vs_{msa_col}"] = {
                        'correlation': float(correlation),
                        'n_samples': len(conf_vals)
                    }
            except Exception as e:
                print(f"Error calculating correlation for {conf_col} vs {msa_col}: {e}")
    
    return correlations

def main():
    parser = argparse.ArgumentParser(
        description="Generate summary statistics for MSA analysis results"
    )
    parser.add_argument(
        "--pmhc", "-p",
        type=Path,
        required=True,
        help="p:MHC metrics parquet file"
    )
    parser.add_argument(
        "--triad", "-t",
        type=Path,
        required=True,
        help="Triad metrics parquet file"
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        required=True,
        help="Output JSON file"
    )
    
    args = parser.parse_args()
    
    # Load data
    print("Loading p:MHC data...")
    pmhc_df = pl.read_parquet(args.pmhc)
    
    print("Loading triad data...")
    triad_df = pl.read_parquet(args.triad)
    
    # Generate summary
    summary = {
        'pmhc_analysis': {
            'n_complexes': len(pmhc_df),
            'summary_statistics': calculate_summary_stats(pmhc_df, 'pmhc'),
            'paired_vs_unpaired_comparison': compare_paired_vs_unpaired(pmhc_df, 'pmhc'),
            'confidence_correlations': analyze_confidence_correlation(pmhc_df)
        },
        'triad_analysis': {
            'n_triads': len(triad_df),
            'summary_statistics': calculate_summary_stats(triad_df, 'triad'),
            'paired_vs_unpaired_comparison': compare_paired_vs_unpaired(triad_df, 'triad'),
            'confidence_correlations': analyze_confidence_correlation(triad_df)
        }
    }
    
    # Save summary
    with open(args.output, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"Summary saved to {args.output}")
    
    # Print key findings
    print("\n=== Key Findings ===")
    print(f"p:MHC complexes analyzed: {summary['pmhc_analysis']['n_complexes']}")
    print(f"Triads analyzed: {summary['triad_analysis']['n_triads']}")
    
    # Print Neff comparisons
    print("\nNeff Comparisons (Paired vs Unpaired):")
    for data_type in ['pmhc_analysis', 'triad_analysis']:
        comparisons = summary[data_type]['paired_vs_unpaired_comparison']
        for key, value in comparisons.items():
            if 'neff_comparison' in key:
                seq_type = key.replace('_neff_comparison', '')
                mean_ratio = value.get('mean_ratio_paired_to_unpaired', 0)
                print(f"  {data_type.split('_')[0].upper()} {seq_type}: Mean paired/unpaired ratio = {mean_ratio:.2f}")

if __name__ == "__main__":
    main()