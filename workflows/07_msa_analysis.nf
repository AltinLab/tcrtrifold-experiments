#!/usr/bin/env nextflow

/*
 * MSA Analysis Pipeline
 * 
 * Analyzes MSA metrics (Neff and interface signal) for p:MHC and triad data
 * to evaluate the influence of paired vs unpaired MSAs on confidence metrics.
 */

nextflow.enable.dsl = 2

/*
 * Define parameters
 */
params.input_pmhc = "${baseDir}/data/test/pmhc/staged/test_pmhc.cleaned.parquet"
params.input_triad = "${baseDir}/data/test/triad/staged/test_triad.cleaned.parquet"
params.msa_dir = "${baseDir}/data/test/msa"
params.output_dir = "${baseDir}/results/msa_analysis"
params.sample_size = null  // Use all data if null
params.help = false

/*
 * Print help message
 */
def helpMessage() {
    log.info """
    MSA Analysis Pipeline
    =====================
    
    This pipeline analyzes MSA metrics for p:MHC and triad data to evaluate
    the influence of paired and unpaired MSAs on confidence metrics.
    
    Usage:
        nextflow run 07_msa_analysis.nf [options]
    
    Options:
        --input_pmhc PATH       Input p:MHC parquet file (default: test data)
        --input_triad PATH      Input triad parquet file (default: test data)
        --msa_dir PATH          Directory containing MSA files (default: test data)
        --output_dir PATH       Output directory (default: results/msa_analysis)
        --sample_size INT       Sample size for testing (optional)
        --help                  Show this help message
    
    Output:
        - pmhc_msa_metrics.parquet: p:MHC data with MSA metrics
        - triad_msa_metrics.parquet: Triad data with MSA metrics
        - msa_analysis_summary.json: Summary statistics
    """
}

if (params.help) {
    helpMessage()
    exit 0
}

/*
 * Log parameters
 */
log.info """
MSA Analysis Pipeline
=====================
Input p:MHC file    : ${params.input_pmhc}
Input triad file    : ${params.input_triad}
MSA directory       : ${params.msa_dir}
Output directory    : ${params.output_dir}
Sample size         : ${params.sample_size ?: 'all data'}
"""

/*
 * Process p:MHC data
 */
process ANALYZE_PMHC_MSA {
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path input_file
    path msa_dir
    val sample_size
    
    output:
    path "pmhc_msa_metrics.parquet", emit: pmhc_metrics
    
    script:
    sample_arg = sample_size ? "--sample-size ${sample_size}" : ""
    """
    python ${baseDir}/workflows/bin/analyze_msa_metrics.py \\
        --input ${input_file} \\
        --msa-dir ${msa_dir} \\
        --output pmhc_msa_metrics.parquet \\
        --data-type pmhc \\
        ${sample_arg}
    """
}

/*
 * Process triad data
 */
process ANALYZE_TRIAD_MSA {
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path input_file
    path msa_dir
    val sample_size
    
    output:
    path "triad_msa_metrics.parquet", emit: triad_metrics
    
    script:
    sample_arg = sample_size ? "--sample-size ${sample_size}" : ""
    """
    python ${baseDir}/workflows/bin/analyze_msa_metrics.py \\
        --input ${input_file} \\
        --msa-dir ${msa_dir} \\
        --output triad_msa_metrics.parquet \\
        --data-type triad \\
        ${sample_arg}
    """
}

/*
 * Generate summary report
 */
process GENERATE_SUMMARY {
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path pmhc_metrics
    path triad_metrics
    
    output:
    path "msa_analysis_summary.json"
    
    script:
    """
    python ${baseDir}/workflows/bin/generate_msa_summary.py \\
        --pmhc ${pmhc_metrics} \\
        --triad ${triad_metrics} \\
        --output msa_analysis_summary.json
    """
}

/*
 * Main workflow
 */
workflow {
    // Create channels for input files
    input_pmhc_ch = Channel.fromPath(params.input_pmhc, checkIfExists: true)
    input_triad_ch = Channel.fromPath(params.input_triad, checkIfExists: true)
    msa_dir_ch = Channel.fromPath(params.msa_dir, type: 'dir', checkIfExists: true)
    
    // Process p:MHC data
    pmhc_metrics = ANALYZE_PMHC_MSA(
        input_pmhc_ch,
        msa_dir_ch,
        params.sample_size
    )
    
    // Process triad data
    triad_metrics = ANALYZE_TRIAD_MSA(
        input_triad_ch,
        msa_dir_ch,
        params.sample_size
    )
    
    // Generate summary
    GENERATE_SUMMARY(
        pmhc_metrics.pmhc_metrics,
        triad_metrics.triad_metrics
    )
}

/*
 * Completion message
 */
workflow.onComplete {
    log.info """
    Pipeline completed!
    Results are available in: ${params.output_dir}
    """
}