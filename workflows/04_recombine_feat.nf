



process JOIN_TRIAD_FEATURES {
    label "tcrtrifold_local"

    publishDir "${params.data_dir}/${params.dset_name}/triad/staged", mode: 'copy'

    input:
    path base_triad_file
    val triad_feature_glob

    output:
    path("*.parquet")

    script:
    """
    #!/usr/bin/env python
    import polars as pl
    from pathlib import Path

    feat_fnames = "${triad_feature_glob.join(',')}".split(',')

    base_triad = pl.read_parquet("${base_triad_file}")
    
    out_triad = base_triad

    for feat_fname in feat_fnames:
        triad_features = pl.read_parquet(feat_fname)
        out_triad = out_triad.join(triad_features, on="job_name")

    out_triad.write_parquet("${base_triad_file.getBaseName()}.feat.parquet")
    """
}

process JOIN_PMHC_FEATURES {
    label "tcrtrifold_local"

    publishDir "${params.data_dir}/${params.dset_name}/pmhc/staged", mode: 'copy'

    input:
    path base_pmhc_file
    val triad_feature_glob

    output:
    path("*.parquet")

    script:
    """
    #!/usr/bin/env python
    import polars as pl
    from pathlib import Path

    feat_fnames = "${triad_feature_glob.join(',')}".split(',')

    base_pmhc = pl.read_parquet("${base_pmhc_file}")
    
    out_pmhc = base_pmhc

    for feat_fname in feat_fnames:
        pmhc_feat = pl.read_parquet(feat_fname)
        out_pmhc = out_pmhc.join(pmhc_feat, on="job_name")

    out_pmhc.write_parquet("${base_pmhc_file.getBaseName()}.feat.parquet")
    """
}

workflow {

    base_triad_channel = Channel.fromPath("${params.data_dir}/${params.dset_name}/triad/features/${params.input_pattern}").collect().map {
        list -> list[0]
    }

    triad_feature_glob = Channel.fromPath("${params.data_dir}/${params.dset_name}/triad/features/${params.input_pattern}").collect()

    JOIN_TRIAD_FEATURES(
        base_triad_channel,
        triad_feature_glob
    )

    base_pmhc_channel = Channel.fromPath("${params.data_dir}/${params.dset_name}/pmhc/features/${params.input_pattern}").collect().map {
        list -> list[0]
    }

    pmhc_feature_glob = Channel.fromPath("${params.data_dir}/${params.dset_name}/pmhc/features/${params.input_pattern}").collect()

    JOIN_PMHC_FEATURES(
        base_pmhc_channel,
        pmhc_feature_glob
    )
}