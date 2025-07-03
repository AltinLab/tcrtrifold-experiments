params.outdir = "$params.data_dir/$params.dset_name/triad"
params.skip_msa = "0"

include { splitParquet } from 'plugin/nf-parquet'
include { INFERENCE_WORKFLOW }from './subworkflows/tgen/af3'
include { SEQ_LIST_TO_FASTA } from './modules/tgen/af3'


workflow {
 
    triad_channel = Channel.fromPath("$params.data_dir/$params.dset_name/triad/staged/*.neg*.parquet").splitParquet()
        .map{
            row -> 
                if (row["mhc_2_seq"] == null) {
                    tuple(
                        [
                            id : row["job_name"],
                            protein_types : ["peptide", "mhc", "tcr", "tcr"],
                        ],
                        [row["peptide"], row["mhc_1_seq"], row["tcr_1_seq"], row["tcr_2_seq"]],
                    )
                }
                else {
                    tuple(
                        [
                            id : row["job_name"],
                            protein_types : ["peptide", "mhc", "mhc", "tcr", "tcr"],
                        ],
                        [row["peptide"], row["mhc_1_seq"], row["mhc_2_seq"], row["tcr_1_seq"], row["tcr_2_seq"]],
                    )
                }
        }

    triad_fasta_channel = SEQ_LIST_TO_FASTA(triad_channel)
    INFERENCE_WORKFLOW(triad_fasta_channel)

}