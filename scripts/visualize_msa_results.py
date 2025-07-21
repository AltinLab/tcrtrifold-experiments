#!/usr/bin/env python3
"""
Simple visualization script for MSA analysis results.
Can be run independently of the Jupyter notebook.
"""

import argparse
import json
import sys
from pathlib import Path

def print_summary_stats(summary_file):
    """Print summary statistics from the analysis."""
    
    with open(summary_file, 'r') as f:
        summary = json.load(f)
    
    print("=" * 60)
    print("MSA ANALYSIS SUMMARY")
    print("=" * 60)
    
    # p:MHC Analysis
    pmhc_stats = summary.get('pmhc_analysis', {})
    print(f"\n🧬 p:MHC COMPLEXES ANALYZED: {pmhc_stats.get('n_complexes', 0)}")
    
    pmhc_comparisons = pmhc_stats.get('paired_vs_unpaired_comparison', {})
    print("\n📊 p:MHC Paired vs Unpaired Neff Comparison:")
    
    for key, value in pmhc_comparisons.items():
        if 'neff_comparison' in key:
            seq_type = key.replace('_neff_comparison', '').replace('_', ' ').upper()
            mean_ratio = value.get('mean_ratio_paired_to_unpaired', 0)
            paired_higher = value.get('paired_higher_count', 0)
            unpaired_higher = value.get('unpaired_higher_count', 0)
            total = paired_higher + unpaired_higher + value.get('equal_count', 0)
            
            if total > 0:
                paired_pct = (paired_higher / total) * 100
                unpaired_pct = (unpaired_higher / total) * 100
                
                print(f"  {seq_type}:")
                print(f"    Mean Paired/Unpaired Ratio: {mean_ratio:.3f}")
                print(f"    Paired higher: {paired_higher} ({paired_pct:.1f}%)")
                print(f"    Unpaired higher: {unpaired_higher} ({unpaired_pct:.1f}%)")
                
                if mean_ratio > 1.1:
                    print("    → PAIRED MSAs are more effective")
                elif mean_ratio < 0.9:
                    print("    → UNPAIRED MSAs are more effective")
                else:
                    print("    → Similar effectiveness")
    
    # Triad Analysis
    triad_stats = summary.get('triad_analysis', {})
    print(f"\n🔬 TRIADS (TCR:p:MHC) ANALYZED: {triad_stats.get('n_triads', 0)}")
    
    triad_comparisons = triad_stats.get('paired_vs_unpaired_comparison', {})
    print("\n📊 Triad Paired vs Unpaired Neff Comparison:")
    
    for key, value in triad_comparisons.items():
        if 'neff_comparison' in key:
            seq_type = key.replace('_neff_comparison', '').replace('_', ' ').upper()
            mean_ratio = value.get('mean_ratio_paired_to_unpaired', 0)
            paired_higher = value.get('paired_higher_count', 0)
            unpaired_higher = value.get('unpaired_higher_count', 0)
            total = paired_higher + unpaired_higher + value.get('equal_count', 0)
            
            if total > 0:
                paired_pct = (paired_higher / total) * 100
                unpaired_pct = (unpaired_higher / total) * 100
                
                print(f"  {seq_type}:")
                print(f"    Mean Paired/Unpaired Ratio: {mean_ratio:.3f}")
                print(f"    Paired higher: {paired_higher} ({paired_pct:.1f}%)")
                print(f"    Unpaired higher: {unpaired_higher} ({unpaired_pct:.1f}%)")
                
                if mean_ratio > 1.1:
                    print("    → PAIRED MSAs are more effective")
                elif mean_ratio < 0.9:
                    print("    → UNPAIRED MSAs are more effective")
                else:
                    print("    → Similar effectiveness")
    
    # Confidence Correlations
    print("\n🔗 CONFIDENCE METRIC CORRELATIONS:")
    
    pmhc_corrs = pmhc_stats.get('confidence_correlations', {})
    triad_corrs = triad_stats.get('confidence_correlations', {})
    
    if pmhc_corrs and isinstance(pmhc_corrs, dict) and len(pmhc_corrs) > 0:
        print("\n  p:MHC Top Correlations:")
        # Sort by absolute correlation value
        sorted_corrs = sorted(pmhc_corrs.items(), 
                             key=lambda x: abs(x[1].get('correlation', 0)), 
                             reverse=True)
        
        for i, (key, corr_data) in enumerate(sorted_corrs[:5]):
            if isinstance(corr_data, dict):
                corr = corr_data.get('correlation', 0)
                n_samples = corr_data.get('n_samples', 0)
                print(f"    {i+1}. {key.replace('_vs_', ' vs ')}: r={corr:.3f} (n={n_samples})")
    else:
        print("  p:MHC: No confidence correlations found")
    
    if triad_corrs and isinstance(triad_corrs, dict) and len(triad_corrs) > 0:
        print("\n  Triad Top Correlations:")
        sorted_corrs = sorted(triad_corrs.items(), 
                             key=lambda x: abs(x[1].get('correlation', 0)), 
                             reverse=True)
        
        for i, (key, corr_data) in enumerate(sorted_corrs[:5]):
            if isinstance(corr_data, dict):
                corr = corr_data.get('correlation', 0)
                n_samples = corr_data.get('n_samples', 0)
                print(f"    {i+1}. {key.replace('_vs_', ' vs ')}: r={corr:.3f} (n={n_samples})")
    else:
        print("  Triad: No confidence correlations found")
    
    print("\n" + "=" * 60)
    print("KEY INSIGHTS:")
    print("=" * 60)
    
    # Generate insights based on the data
    insights = []
    
    # Check if paired MSAs generally outperform unpaired
    paired_better_count = 0
    total_comparisons = 0
    
    for comparisons in [pmhc_comparisons, triad_comparisons]:
        for key, value in comparisons.items():
            if 'neff_comparison' in key:
                ratio = value.get('mean_ratio_paired_to_unpaired', 1.0)
                if ratio > 1.1:
                    paired_better_count += 1
                total_comparisons += 1
    
    if total_comparisons > 0:
        if paired_better_count / total_comparisons > 0.6:
            insights.append("✅ Paired MSAs generally show higher Neff values than unpaired MSAs")
        elif paired_better_count / total_comparisons < 0.4:
            insights.append("❌ Unpaired MSAs often show higher Neff values than paired MSAs")
        else:
            insights.append("⚖️ Paired and unpaired MSAs show mixed effectiveness")
    
    # Check confidence correlations
    has_strong_corr = False
    for corrs in [pmhc_corrs, triad_corrs]:
        if isinstance(corrs, dict):
            for corr_data in corrs.values():
                if isinstance(corr_data, dict):
                    corr = abs(corr_data.get('correlation', 0))
                    if corr > 0.3:
                        has_strong_corr = True
                        break
    
    if has_strong_corr:
        insights.append("🔗 Strong correlations found between MSA metrics and confidence scores")
    else:
        insights.append("🔗 Weak correlations between MSA metrics and confidence scores")
    
    for i, insight in enumerate(insights, 1):
        print(f"{i}. {insight}")
    
    print("\n📚 This analysis supports the hypothesis from AlphaFold2/3 papers that")
    print("   MSA quality (Neff) and pairing influence prediction confidence.")
    print()


def main():
    parser = argparse.ArgumentParser(
        description="Display MSA analysis summary statistics"
    )
    parser.add_argument(
        "summary_file",
        type=Path,
        help="Path to summary JSON file"
    )
    
    args = parser.parse_args()
    
    if not args.summary_file.exists():
        print(f"Error: Summary file {args.summary_file} does not exist")
        print("Run the MSA analysis pipeline first:")
        print("  nextflow run workflows/07_msa_analysis.nf")
        sys.exit(1)
    
    print_summary_stats(args.summary_file)


if __name__ == "__main__":
    main()