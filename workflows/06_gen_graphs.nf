

process GEN_GRAPHS {
    label "tcrtrifold_local"

    publishDir "${params.data_dir}/${params.dset_name}/triad/graphs", mode: 'copy'

    input:
    path input_pq
    path inference_dir

    output:
    path("*")

    script:
    """
    gen_graphs.py \\
        --input_parquet ${input_pq} \\
        --inference_dir ${inference_dir} \\
        --inference_type ${params.inf_type} \\
        --output_path "${input_pq.getSimpleName()}_${params.inf_type}"
    """

}

workflow {
    clean_pq = Channel.fromPath("${params.data_dir}/${params.dset_name}/triad/staged/*.neg*.parquet")
    
    if (params.inf_type == "af3") {
        inf_dir = Channel.fromPath("${params.data_dir}/${params.dset_name}/triad/inference")
    }
    else {
        inf_dir = Channel.fromPath("${params.data_dir}/${params.dset_name}/triad/predictions")
    }
    

    GEN_GRAPHS(clean_pq, inf_dir)

}