params.outdir = "$params.data_dir/$params.dset_name/triad"

include { splitParquet } from 'plugin/nf-parquet'
include { SEQ_LIST_TO_FASTA } from './modules/tgen/af3'
include { FASTA_TO_YAML; BOLTZ_INFERENCE } from './modules/local/boltz'


workflow {
 
    triad_channel = Channel.fromPath("$params.data_dir/$params.dset_name/triad/staged/*.neg*.parquet").splitParquet()
        .map{
            row -> 
                if (row["mhc_2_seq"] == null) {
                    tuple(
                        [
                            id : row["job_name"],
                            segids : ["A", "B", "D", "E"]
                        ],
                        [row["peptide"], row["mhc_1_seq"], row["tcr_1_seq"], row["tcr_2_seq"]],
                    )
                }
                else {
                    tuple(
                        [
                            id : row["job_name"],
                            segids : ["A", "B", "C", "D", "E"]
                        ],
                        [row["peptide"], row["mhc_1_seq"], row["mhc_2_seq"], row["tcr_1_seq"], row["tcr_2_seq"]],
                    )
                }
        }

    triad_fasta_channel = SEQ_LIST_TO_FASTA(triad_channel)
    triad_yaml_channel = FASTA_TO_YAML(triad_fasta_channel)
    BOLTZ_INFERENCE(triad_yaml_channel)

}