
process GEN_NEGATIVES {
  label "tcrtrifold_heavy"

  publishDir "${params.data_dir}/${params.dset_name}/triad/staged", mode: 'copy'

  input:
  path base_df
  val supp_dfs

  output:
  path("*neg*.parquet")
  path("*discard*.parquet"), optional: true

  script:
  """
  gen_negatives.py \\
    --base_df ${base_df} \\
    --supp_dfs ${supp_dfs.join(",")} \\
    --neg_depth ${params.neg_depth} \\
    -d "${base_df.getSimpleName()}.discard.parquet" \\
    -o "${base_df.getSimpleName()}.neg.parquet"
  """
}


workflow {

    // groovy list of dset names
    cleaned_neg_dset_names = params.negs_from.split(',').collect { neg_dset_name -> neg_dset_name.replace(' ', '')}
    println(cleaned_neg_dset_names)

    per_neg_dset_ch = cleaned_neg_dset_names.collect { dset ->
        Channel.fromPath("${params.data_dir}/${dset}/triad/staged/*cleaned*.parquet")
    }

    if (per_neg_dset_ch.size() > 1) {
        all_neg_dset_ch = per_neg_dset_ch[0].concat(*per_neg_dset_ch[1..(per_neg_dset_ch.size() - 1)]).toList()
    }
    else {
        all_neg_dset_ch = per_neg_dset_ch[0].toList()
    }

    // cleaned_neg_dset_channels = []

    // for (neg_dset_name in cleaned_neg_dset_names) {
    //   cleaned_neg_dset_channels.add(Channel.fromPath("$params.data_dir/$neg_dset_name/triad/staged/*cleaned*.parquet"))
    // }

    // cleaned_neg_dset_channels = 
    //                 .flatMap { dset_name ->                                                  
    //                 Channel.fromPath("${params.data_dir}/${dset_name}/triad/staged/*cleaned*.parquet")
    //                 }.collect()

    // negs_from_list = Channel
    //   .from(params.negs_from.split(',') as List)
    //   .flatMap { dset ->
    //     Channel.fromPath("$params.data_dir/$dset/triad/staged/*cleaned*.parquet")
    //   }
    //   .toList()

    GEN_NEGATIVES(Channel.fromPath("${params.data_dir}/${params.dset_name}/triad/staged/*cleaned*.parquet"),
                    all_neg_dset_ch)

}