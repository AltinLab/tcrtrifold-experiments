#!/usr/bin/env nextflow

/*
 * Pipeline for calculating interface similarity between complexes and their MSA alignments
 */

nextflow.enable.dsl=2

params.data_dir = "data"
params.dset_name = "test"
params.outdir = "$params.data_dir/$params.dset_name"

workflow {
    
    // Process p:MHC complexes
    pmhc_channel = Channel.fromPath("$params.data_dir/$params.dset_name/pmhc/staged/*.cleaned*.parquet")
        .map { file -> tuple("pmhc", file) }
    
    // Process triad complexes
    triad_channel = Channel.fromPath("$params.data_dir/$params.dset_name/triad/staged/*.cleaned*.parquet")
        .map { file -> tuple("triad", file) }
    
    // Combine channels
    all_complexes = pmhc_channel.concat(triad_channel)
    
    // Calculate interface similarities
    CALCULATE_INTERFACE_SIMILARITY(all_complexes)
}

process CALCULATE_INTERFACE_SIMILARITY {
    publishDir "$params.outdir/interface_similarity", mode: 'copy'
    label 'tcrtrifold_local'
    
    input:
    tuple val(complex_type), path(parquet_file)
    
    output:
    path("*_interface_similarity.parquet")
    
    script:
    """
    #!/usr/bin/env python3
    
    import polars as pl
    from pathlib import Path
    from tcrtrifold.interface_similarity import add_interface_similarity_features
    
    # Read the parquet file
    df = pl.read_parquet("${parquet_file}")
    
    # Set up MSA directory
    msa_dir = Path("${params.data_dir}/${params.dset_name}/msa")
    
    # Add interface similarity features
    df_with_features = add_interface_similarity_features(
        df, 
        msa_dir, 
        complex_type="${complex_type}"
    )
    
    # Save output
    output_file = "${parquet_file.baseName}_interface_similarity.parquet"
    df_with_features.write_parquet(output_file)
    
    print(f"Processed {df.shape[0]} rows, added {len(df_with_features.columns) - len(df.columns)} interface similarity features")
    """
}