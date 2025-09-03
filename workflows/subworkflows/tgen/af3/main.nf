/*
Parameters:

MSA:
- force_update_msa: If true, force update the MSA even if it already exists

Inference:
- compress_inf: If true, compress the inference directory after inference using mdaf3
- seeds: string of comma-separated integers to use as random seeds for inference
- collate_inf_size: Batch size for collating inference jobs
- inference_dir: Directory to store inference results. Defaults to "${params.outdir}/inference"
- check_inf_exists: If true, check if inference directories already exist (in inference_dir) before running inference for a job
- save_embeddings: If true, save embeddings during inference
*/




include { FILT_FORMAT_MSA;
            RUN_MSA;
            COMPOSE_INFERENCE_JSON;
            BATCHED_INFERENCE;
            INFERENCE;
            CLEAN_INFERENCE_DIR} from '../../../modules/tgen/af3'





workflow MSA_WORKFLOW {
    /*
    Arguments:
    - meta_fasta: Channel of tuples containing metadata and fasta files. Metadata should contain:
        - id: Unique identifier for the sequence
        - protein_types: List of protein types for each chain in the fasta (one of ["peptide", "mhc", "tcr", "any"])
        - segids (optional): List of segment IDs to label chains in the output structure
        - skip_msa (optional): List of integers to skip MSA for specific indices in the input fasta file

    Returns:
    - new_msa_list: List of newly completed MSAs. Can be used as a single token representing MSA completion
    */

    take:
    meta_fasta

    main:
    FILT_FORMAT_MSA(meta_fasta)
    RUN_MSA(FILT_FORMAT_MSA.out)
    // STORE_MSA(RUN_MSA.out)

    emit:
    new_msa_list = RUN_MSA.out.toList()

}



workflow INFERENCE_WORKFLOW {
    /*
    Arguments:
    - meta_fasta: Channel of tuples containing metadata and fasta files. Metadata should contain:
        - id: Unique identifier for the sequence
        - protein_types: List of protein types (one of ["peptide", "mhc", "tcr", "any"])
        - segids (optional): List of segment IDs to label chains in the output structure
    - inf_dir: Directory to check for existing inference results

    Returns:
    - new_inf_list: List of completed inferences. Can be used as a single token representing inference completion
    */

    take:
    meta_fasta
    inf_dir

    main:

    json = COMPOSE_INFERENCE_JSON(meta_fasta, inf_dir)

    batched_json = json.collate(params.collate_inf_size).map { batch ->
        def allMeta    = batch.collect { it[0] }
        def allSeqLists = batch.collect { it[1] }
        tuple(allMeta, allSeqLists)
    }

    inference = BATCHED_INFERENCE(batched_json)

    inference = inference.map { meta, inf_dirs ->
        def listOut = (inf_dirs instanceof List) ? inf_dirs : [ inf_dirs ]
        tuple(meta, listOut)
    }

    inference = inference.flatMap { metas, inf_dirs ->
        metas.indices.collect { idx ->
            tuple( metas[idx], inf_dirs[idx] )
        }
    }

    if (params.compress_inf == true) {
        // Clean up inference directory
        inference = CLEAN_INFERENCE_DIR(inference)
    }


    emit:
    new_meta_inf = inference

}

workflow UNBATCHED_INFERENCE_WORKFLOW {
    /*
    Arguments:
    - meta_fasta: Channel of tuples containing metadata and fasta files. Metadata should contain:
        - id: Unique identifier for the sequence
        - protein_types: List of protein types (one of ["peptide", "mhc", "tcr", "any"])
        - segids (optional): List of segment IDs to label chains in the output structure
    - inf_dir: Directory to check for existing inference results

    Returns:
    - new_inf_list: List of completed inferences. Can be used as a single token representing inference completion
    */

    take:
    meta_fasta
    inf_dir

    main:

    json = COMPOSE_INFERENCE_JSON(meta_fasta, inf_dir)

    inference = INFERENCE(json)

    if (params.compress_inf == true) {
        // Clean up inference directory
        inference = CLEAN_INFERENCE_DIR(inference)
    }


    emit:
    new_meta_inf = inference

}


