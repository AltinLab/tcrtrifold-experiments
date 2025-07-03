include { splitParquet } from 'plugin/nf-parquet'
include { MSA_WORKFLOW } from './subworkflows/tgen/af3'
include { SEQ_LIST_TO_FASTA } from './modules/tgen/af3'

workflow {

    mhc_1_channel = Channel.fromPath("$params.data_dir/$params.dset_name/triad/staged/*.neg*.parquet").splitParquet()
    .map{
        row -> 
            tuple(
                [
                    id : "mhc_1",
                    protein_type : "mhc",
                ],
                [row["mhc_1_seq"]],
            )
    }.unique()

    mhc_2_channel = Channel.fromPath("$params.data_dir/$params.dset_name/triad/staged/*.neg*.parquet").splitParquet()
    .filter { row -> 
                row["mhc_2_seq"] != null
            }
    .map{
        row -> 
            tuple(
                [
                    id : "mhc_2",
                    protein_type : "mhc",
                ],
                [row["mhc_2_seq"]],
            )
    }.unique()

    tcr_1_channel = Channel.fromPath("$params.data_dir/$params.dset_name/triad/staged/*.neg*.parquet").splitParquet()
    .map{
        row -> 
            tuple(
                [
                    id : "tcr_1",
                    protein_type : "tcr",
                ],
                [row["tcr_1_seq"]],
            )
    }.unique()

    tcr_2_channel = Channel.fromPath("$params.data_dir/$params.dset_name/triad/staged/*.neg*.parquet").splitParquet()
    .map{
        row -> 
            tuple(
                [
                    id : "tcr_2",
                    protein_type : "tcr",
                ],
                [row["tcr_2_seq"]],
            )
    }.unique()

    all_proteins = mhc_1_channel.concat(mhc_2_channel, tcr_1_channel, tcr_2_channel)

    all_proteins_fasta_channel = SEQ_LIST_TO_FASTA(all_proteins)
 
    MSA_WORKFLOW(all_proteins_fasta_channel)
}