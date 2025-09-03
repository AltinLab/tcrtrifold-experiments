include { SEQ_LIST_TO_FASTA } from '../../../modules/tgen/af3'
include { MSA_WORKFLOW; INFERENCE_WORKFLOW; UNBATCHED_INFERENCE_WORKFLOW } from '../../tgen/af3'
include { splitParquet } from 'plugin/nf-parquet'

def hash_from(String seq) {
    def md = java.security.MessageDigest.getInstance("SHA-256")
    md.update(seq.getBytes('UTF-8'))      // use getBytes(...) rather than .bytes
    def bytes = md.digest()
    bytes.collect { String.format('%02x', it) }.join()
}

process NOOP_DEP {

    label "process_local"

    input:
    path input_val
    val msa_done_key

    output:
    path "*", includeInputs: true

    script:
    """
    true
    """
}

workflow MSA_FROM_TRIAD_PARQUET {
    take:
    triad_parquet

    main:

    mhc_1_channel = triad_parquet.splitParquet()
    .map{
        row -> 

        def sha256 = hash_from(row["mhc_1_seq"])


        tuple(
            [
                id : sha256,
                protein_type : "mhc",
            ],
            [row["mhc_1_seq"]],
        )
    }.unique()

    mhc_2_channel = triad_parquet.splitParquet()
    .filter { row -> 
                row["mhc_2_seq"] != null
            }
    .map{
        row -> 

        def sha256 = hash_from(row["mhc_2_seq"])
            tuple(
                [
                    id : sha256,
                    protein_type : "mhc",
                ],
                [row["mhc_2_seq"]],
            )
    }.unique()

    tcr_1_channel = triad_parquet.splitParquet()
    .map{
        row -> 

        def sha256 = hash_from(row["tcr_1_seq"])
        tuple(
            [
                id : sha256,
                protein_type : "tcr",
            ],
            [row["tcr_1_seq"]],
        )
    }.unique()

    tcr_2_channel = triad_parquet.splitParquet()
    .map{
        row -> 

        def sha256 = hash_from(row["tcr_2_seq"])
        tuple(
            [
                id : sha256,
                protein_type : "tcr",
            ],
            [row["tcr_2_seq"]],
        )
    }.unique()

    all_proteins = mhc_1_channel.concat(mhc_2_channel, tcr_1_channel, tcr_2_channel)

    all_proteins_fasta_channel = SEQ_LIST_TO_FASTA(all_proteins)
 
    MSA_WORKFLOW(all_proteins_fasta_channel)


    emit:
    new_msa_list = MSA_WORKFLOW.out.new_msa_list

}

workflow INFERENCE_FROM_TRIAD_PARQUET {

    take:
    triad_parquet
    triad_inf_dir
    msa_done_key

    main:

    triad_parquet = NOOP_DEP(triad_parquet, msa_done_key)
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
    INFERENCE_WORKFLOW(triad_fasta_channel, triad_inf_dir)

    emit:
    new_meta_inf = INFERENCE_WORKFLOW.out.new_meta_inf
}


workflow UNBATCHED_INFERENCE_FROM_TRIAD_PARQUET {

    take:
    triad_parquet
    triad_inf_dir
    msa_done_key

    main:

    triad_parquet = NOOP_DEP(triad_parquet, msa_done_key)
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
    UNBATCHED_INFERENCE_WORKFLOW(triad_fasta_channel, triad_inf_dir)

    emit:
    new_meta_inf = UNBATCHED_INFERENCE_WORKFLOW.out.new_meta_inf
}


// workflow MSA_FROM_PMHC_PARQUET {
//     take:
//     triad_parquet

//     main:

//     mhc_1_channel = triad_parquet.splitParquet()
//     .map{
//         row -> 

//             def sha256 = hash_from(row["mhc_1_seq"])
//             tuple(
//                 [
//                     id : sha256,
//                     protein_type : "mhc",
//                 ],
//                 [row["mhc_1_seq"]],
//             )
//     }.unique()

//     mhc_2_channel = triad_parquet.splitParquet()
//     .filter { row -> 
//                 row["mhc_2_seq"] != null
//             }
//     .map{

        
//         row -> 

//             def sha256 = hash_from(row["mhc_2_seq"])
//             tuple(
//                 [
//                     id : sha256,
//                     protein_type : "mhc",
//                 ],
//                 [row["mhc_2_seq"]],
//             )
//     }.unique()

//     tcr_1_channel = triad_parquet.splitParquet()
//     .map{
//         row -> 

//             def sha256 = hash_from(row["tcr_1_seq"])
//             tuple(
//                 [
//                     id : sha256,
//                     protein_type : "tcr",
//                 ],
//                 [row["tcr_1_seq"]],
//             )
//     }.unique()

//     tcr_2_channel = triad_parquet.splitParquet()
//     .map{
//         row -> 

//             def sha256 = hash_from(row["tcr_2_seq"])
//             tuple(
//                 [
//                     id : sha256,
//                     protein_type : "tcr",
//                 ],
//                 [row["tcr_2_seq"]],
//             )
//     }.unique()

//     all_proteins = mhc_1_channel.concat(mhc_2_channel, tcr_1_channel, tcr_2_channel)

//     all_proteins_fasta_channel = SEQ_LIST_TO_FASTA(all_proteins)
 
//     MSA_WORKFLOW(all_proteins_fasta_channel)


//     emit:
//     new_msa_list = MSA_WORKFLOW.out.new_msa_list

// }

workflow INFERENCE_FROM_PMHC_PARQUET {

    take:
    pmhc_parquet
    pmhc_inf_dir
    msa_done_key

    main:

    pmhc_parquet = NOOP_DEP(pmhc_parquet, msa_done_key)
    pmhc_channel = pmhc_parquet.splitParquet()
    .map{
        row -> 
            if (row["mhc_2_seq"] == null) {
                tuple(
                    [
                        id : row["job_name"],
                        protein_types : ["peptide", "mhc"],
                        segids : ["A", "B"],
                        skip_msa : [0]
                    ],
                    [row["peptide"], row["mhc_1_seq"]],
                )
            }
            else {
                tuple(
                    [
                        id : row["job_name"],
                        protein_types : ["peptide", "mhc", "mhc"],
                        segids : ["A", "B", "C"],
                        skip_msa : [0]
                    ],
                    [row["peptide"], row["mhc_1_seq"], row["mhc_2_seq"]],
                )
            }
    }

    pmhc_fasta_channel = SEQ_LIST_TO_FASTA(pmhc_channel)
    INFERENCE_WORKFLOW(pmhc_fasta_channel, pmhc_inf_dir)

    emit:
    new_meta_inf = INFERENCE_WORKFLOW.out.new_meta_inf
}


workflow UNBATCHED_INFERENCE_FROM_PMHC_PARQUET {

    take:
    pmhc_parquet
    pmhc_inf_dir
    msa_done_key

    main:

    pmhc_parquet = NOOP_DEP(pmhc_parquet, msa_done_key)
    pmhc_channel = pmhc_parquet.splitParquet()
    .map{
        row -> 
            if (row["mhc_2_seq"] == null) {
                tuple(
                    [
                        id : row["job_name"],
                        protein_types : ["peptide", "mhc"],
                        segids : ["A", "B"],
                        skip_msa : [0]
                    ],
                    [row["peptide"], row["mhc_1_seq"]],
                )
            }
            else {
                tuple(
                    [
                        id : row["job_name"],
                        protein_types : ["peptide", "mhc", "mhc"],
                        segids : ["A", "B", "C"],
                        skip_msa : [0]
                    ],
                    [row["peptide"], row["mhc_1_seq"], row["mhc_2_seq"]],
                )
            }
    }

    pmhc_fasta_channel = SEQ_LIST_TO_FASTA(pmhc_channel)
    UNBATCHED_INFERENCE_WORKFLOW(pmhc_fasta_channel, pmhc_inf_dir)

    emit:
    new_meta_inf = UNBATCHED_INFERENCE_WORKFLOW.out.new_meta_inf
}

