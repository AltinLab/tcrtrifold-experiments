process FASTA_TO_YAML {
    label "boltz_local"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.yaml"), optional: true

    script:
    def segids = (meta.containsKey('segids')) ? "--segids ${meta.segids.join(',')}" : ''
    def skip_msa_arg = (params.skip_msa != null) ? "--skip_msa ${params.skip_msa}" : ''
    def check_inf_exists = params.check_inf_exists ? """
    if [ -d "${params.outdir}/predictions/${meta.id}" ]; then
        echo "Skipping ${meta.id}"
        exit 0
    fi
    """ : ''
    """
    $check_inf_exists

    compose_inference_YAML.py \\
        -jn "${meta.id}" \\
        --fasta_path ${fasta} \\
        ${skip_msa_arg} \\
        ${segids}
    """
}

process BOLTZ_INFERENCE {
    label "boltz_gpu"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(meta), path(yaml)

    output:
    tuple val(meta), path("predictions/*")


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

    mv boltz_results_${yaml.getSimpleName()}/predictions .
    """
}