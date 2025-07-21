params.from_true_struct = false

process TCRDOCK_GEOM {
    label "tcrdock"
    publishDir "${params.data_dir}/${params.dset_name}/triad/staged", mode: 'copy'

    input:
    path input_pq
    path top_dir

    output:
    path("*.parquet")

    script:
    def from_true_struct_arg = params.from_true_struct ? "--from_true_struct" : ""
    def output_fname =  params.from_true_struct ? "${input_pq.getSimpleName()}.true_tcrdock.parquet" : "${input_pq.getSimpleName()}.${params.inf_type}_tcrdock.parquet"
    def inf_type_arg = params.from_true_struct ? "" : "--inference_type ${params.inf_type}"
    """
    tcrdock_geom.py \
        --input_parquet ${input_pq} \\
        --topology_path ${top_dir} \\
        ${from_true_struct_arg} \\
        ${inf_type_arg} \\
        -o ${output_fname}
    """
}


workflow {
    
    if (params.from_true_struct) {
        input_parquet = Channel.fromPath("$params.data_dir/$params.dset_name/triad/staged/*cleaned*.parquet")
        topology_path = Channel.fromPath("$params.data_dir/$params.dset_name/triad/cleaned_pdb")
    }
    else {
        input_parquet = Channel.fromPath("$params.data_dir/$params.dset_name/triad/staged/*neg*.parquet")

        if (params.inf_type == "af3") {
            topology_path = Channel.fromPath("$params.data_dir/$params.dset_name/triad/inference")
        }
        else {
            topology_path = Channel.fromPath("$params.data_dir/$params.dset_name/triad/predictions")
        }
        
    }

    TCRDOCK_GEOM(
        input_parquet,
        topology_path
    )
}