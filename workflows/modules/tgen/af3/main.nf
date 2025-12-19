process SEQ_LIST_TO_FASTA {
    label "process_local"

    input:
      tuple val(meta), val(seq_list)

    output:
      tuple val(meta), path("af3-single-chain-msa.fasta")

    script:
    """
    filename="af3-single-chain-msa.fasta"
    : > "\$filename"

    i=1
    for seq in ${seq_list.join(' ')}; do
        echo ">\$i"   >> "\$filename"
        echo "\$seq"  >> "\$filename"
        i=\$(( i + 1 ))
    done
    """
}

// process FILTER_MISSING_MSA {
//     label "process_local"
//     tag "${meta.protein_type}-${meta.id}"

//     input:
//     tuple val(meta), path(fasta)

//     output:
//     tuple val(meta), path("*.fasta"), optional: true


//     script:
//     """
//     module load singularity

//     fname=\$(uuidgen).csv

//     export SINGULARITYENV_VAST_S3_ACCESS_KEY_ID="\$VAST_S3_ACCESS_KEY_ID"
//     export SINGULARITYENV_VAST_S3_SECRET_ACCESS_KEY="\$VAST_S3_SECRET_ACCESS_KEY"

//     singularity exec --nv \\
//         -B /home,/scratch,/tgen_labs --cleanenv \\
//         /tgen_labs/altin/alphafold3/containers/msa-db.sif \\
//         python ${moduleDir}/resources/usr/bin/filter_missing_msa.py \\
//             -t "${meta.protein_type}" \\
//             -f "$fasta" \\
//             -o "${fasta.getSimpleName()}.filt.json"
//     """
// }

// process COMPOSE_EMPTY_MSA_JSON {
//     label "process_local"
//     tag "${meta.protein_type}-${meta.id}"

//     input:
//     tuple val(meta), path(fasta)

//     output:
//     tuple val(meta), path("*.json")

//     script:
//     """
//     generate_single_JSON.py \\
//         --fasta "$fasta" \\
//         --job_name "${meta.id}"
//     """
//  }

process FILT_FORMAT_MSA {
    label "tcrtrifold_local"
    tag "${meta.protein_type}-${meta.id}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.json"), optional: true


    script:
    def force = params.force_update_msa ? "--force" : ''
    """
   filt_format_msa.py \\
        --protein_type "${meta.protein_type}" \\
        --job_name "${meta.id}" \\
        --fasta "$fasta" \\
        --msa_cache_dir "${params.msa_cache_dir}" \\
        ${force} 
    """
}

process RUN_MSA {
    label "alphafold3_msa"

    memory { "${ Math.min(512, 64 * Math.pow(2, task.attempt - 1)) }GB" }
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 5

    storeDir "${params.msa_cache_dir}/${meta.protein_type}"

    tag "${meta.protein_type}-${meta.id}"

    input:
    tuple val(meta), path(json)

    output:
    tuple val(meta), path("${meta.id}.json")

    script:
    """
    python /app/alphafold/run_alphafold.py \\
        --json_path=$json \\
        --model_dir=${params.af3_model_dir} \\
        --db_dir=${params.af3_db_dir} \\
        --output_dir=. \\
        --norun_inference

    mv ${meta.id}/${meta.id}_data.json ${meta.id}.json
    """
}



// process STORE_MSA {
//     label "process_local"
//     errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
//     maxRetries 5
//     tag "${meta.protein_type}-${meta.id}"
    
//     input:
//     tuple val(meta), path(json)

//     output:
//     tuple val(meta), path(json)

//     script:
//     """
//     module load singularity

//     export SINGULARITYENV_VAST_S3_ACCESS_KEY_ID="\$VAST_S3_ACCESS_KEY_ID"
//     export SINGULARITYENV_VAST_S3_SECRET_ACCESS_KEY="\$VAST_S3_SECRET_ACCESS_KEY"
    
//     singularity exec --nv \\
//         -B /home,/scratch,/tgen_labs --cleanenv \\
//         /tgen_labs/altin/alphafold3/containers/msa-db.sif \\
//         python ${moduleDir}/resources/usr/bin/store_msa.py \\
//             -t "${meta.protein_type}" \\
//             -j "$json"
//     """
// }



process COMPOSE_INFERENCE_JSON {
    label "tcrtrifold_local"
    // errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    // maxRetries 5
    tag "${meta.id}"

    input:
    tuple val(meta), path(fasta)
    // val because we dont want this to be resolved to relative
    val(inf_dir)

    output:
    tuple val(meta), path("*.json"), optional: true


    script:
    def seeds = params.seeds ? "--seeds ${params.seeds}" : ''
    def segids = (meta.containsKey('segids')) ? "--segids ${meta.segids.join(',')}" : ''
    def check_inf_exists = params.check_inf_exists ? "--check_inf_exists" : ''
    def skip_msa_arg = meta.containsKey("skip_msa") ? "--skip_msa ${meta.skip_msa.join(',')}" : ''
    """
    
    
    compose_inference_JSON.py \\
        --job_name "${meta.id}" \\
        --fasta "$fasta" \\
        --protein_types "${meta.protein_types.join(',')}" \\
        --msa_cache_dir "${params.msa_cache_dir}" \\
        --inf_dir "$inf_dir" \\
        ${segids} \\
        ${skip_msa_arg} \\
        ${seeds} \\
        ${check_inf_exists}
    """
 }

process BATCHED_INFERENCE {
    tag "batched_inference"
    label "alphafold3_inference"
    // if (params.compress_inf == false) {
    //     publishDir "${params.outdir}", mode: 'copy'
    // }

    input:
    tuple val(batched_meta), path(batched_json)

    output:
    tuple val(batched_meta), path("*")

    script:
    def save_embeddings = params.save_embeddings ? "--save_embeddings" : ''
    """
    python /app/alphafold/run_alphafold.py \\
        --input_dir=. \\
        --model_dir=$params.af3_model_dir \\
        --db_dir=$params.af3_db_dir \\
        --output_dir=. \\
        --norun_data_pipeline \\
        ${save_embeddings} \\
        --num_diffusion_samples=1
    """
}


process INFERENCE {
    tag "inference"
    label "alphafold3_inference"

    input:
    tuple val(meta), path(json)

    output:
    tuple val(meta), path("*")

    script:
    def save_embeddings = params.save_embeddings ? "--save_embeddings" : ''
    """
    python /app/alphafold/run_alphafold.py \\
        --json_path=$json \\
        --model_dir=$params.af3_model_dir \\
        --db_dir=$params.af3_db_dir \\
        --output_dir=. \\
        --norun_data_pipeline \\
        ${save_embeddings} \\
        --num_diffusion_samples=1
    """
}

process CLEAN_INFERENCE_DIR {
    label "tcrtrifold_local"
    tag "clean_inference"
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 5
    
    // publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(meta), path(inference_dir)

    output:
    tuple val(meta), path("*", includeInputs: true)

    script:
    """
    clean_inference_dir.py \\
        -i $inference_dir
    """
}
