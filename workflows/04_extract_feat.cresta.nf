params.outdir = "$params.data_dir/$params.dset_name"

include { EXTRACT_TRIAD_CONF_FEAT } from './modules/local/extract_feat'


workflow {

    neg_pq = Channel.fromPath("${params.data_dir}/${params.dset_name}/triad/staged/*.neg*.parquet")
    triad_inf_dir = Channel.fromPath("${params.data_dir}/${params.dset_name}/triad/inference")

    EXTRACT_TRIAD_CONF_FEAT(neg_pq, triad_inf_dir)

}