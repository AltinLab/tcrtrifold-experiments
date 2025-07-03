params.outdir = "$params.data_dir/$params.dset_name/pmhc"

include { splitParquet } from 'plugin/nf-parquet'
include { INFERENCE_WORKFLOW }from './subworkflows/tgen/af3'
include { SEQ_LIST_TO_FASTA } from './modules/tgen/af3'


workflow {
 
    pmhc_channel = Channel.fromPath("$params.data_dir/$params.dset_name/pmhc/staged/*.cleaned*.parquet").splitParquet()
        .map{
            row -> 
                if (row["mhc_2_seq"] == null) {
                    tuple(
                        [
                            id : row["job_name"],
                            protein_types : ["peptide", "mhc"],
                        ],
                        [row["peptide"], row["mhc_1_seq"]],
                    )
                }
                else {
                    tuple(
                        [
                            id : row["job_name"],
                            protein_types : ["peptide", "mhc", "mhc"],
                        ],
                        [row["peptide"], row["mhc_1_seq"], row["mhc_2_seq"]],
                    )
                }
        }

    pmhc_fasta_channel = SEQ_LIST_TO_FASTA(pmhc_channel)
    INFERENCE_WORKFLOW(pmhc_fasta_channel)

}