
params.outdir = "${params.data_dir}/iedb_I"

include { MEAN_TCR_PMHC_PAE } from './modules/local/extract_feat/main.nf'

workflow {

    triad = Channel.fromPath("${params.data_dir}/iedb_I/triad/staged/iedb_I.triad.thresh_3.10x_neg.base.parquet")

    MEAN_TCR_PMHC_PAE(
        triad
    )

}