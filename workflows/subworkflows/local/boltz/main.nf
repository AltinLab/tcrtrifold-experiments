include { SEQ_LIST_TO_FASTA } from '../../../modules/tgen/af3'

process FASTA_TO_YAML {
    label "boltz_local"

    input:
    tuple val(meta), path(fasta)
    val inf_dir

    output:
    tuple val(meta), path("*.yaml"), optional: true

    script:
    def segids = (meta.containsKey('segids')) ? "--segids ${meta.segids.join(',')}" : ''
    def skip_msa_arg = meta.containsKey("skip_msa") ? "--skip_msa ${meta.skip_msa.join(',')}" : ''
    """
    compose_inference_YAML.py \\
        -jn "${meta.id}" \\
        --fasta_path ${fasta} \\
        --inf_dir "$inf_dir" \\
        ${skip_msa_arg} \\
        ${segids}
    """
}

process BOLTZ_INFERENCE {
    label "boltz_gpu"

    input:
    tuple val(meta), path(yaml)

    output:
    tuple val(meta), path("${meta.id}")


    script:
    """
    boltz predict \\
        ${yaml} \\
        --cache ${params.boltz_cache} \\
        --use_msa_server \\
        --diffusion_samples 5 \\
        --recycling_steps 5 \\
        --write_full_pae \\
        --write_full_pde

    mv boltz_results_${meta.id}/predictions/${meta.id} .
    """
}

workflow BOLTZ_FROM_TRIAD_PARQUET {

    take:
    triad_parquet
    triad_inf_dir

    main:

    triad_channel = triad_parquet.splitParquet()
    .map{
        row -> 
            if (row["mhc_2_seq"] == null) {
                tuple(
                    [
                        id : row["job_name"],
                        protein_types : ["peptide", "mhc", "tcr", "tcr"],
                        segids : ["A", "B", "D", "E"],
                        skip_msa : [0]
                    ],
                    [row["peptide"], row["mhc_1_seq"], row["tcr_1_seq"], row["tcr_2_seq"]],
                )
            }
            else {
                tuple(
                    [
                        id : row["job_name"],
                        protein_types : ["peptide", "mhc", "mhc", "tcr", "tcr"],
                        segids : ["A", "B", "C", "D", "E"],
                        skip_msa : [0]
                    ],
                    [row["peptide"], row["mhc_1_seq"], row["mhc_2_seq"], row["tcr_1_seq"], row["tcr_2_seq"]],
                )
            }
    }

    triad_fasta_channel = SEQ_LIST_TO_FASTA(triad_channel)
    triad_yaml_channel = FASTA_TO_YAML(triad_fasta_channel, triad_inf_dir)
    
    traid_meta_inf = BOLTZ_INFERENCE(triad_yaml_channel)

    emit:
    new_meta_inf = traid_meta_inf
}