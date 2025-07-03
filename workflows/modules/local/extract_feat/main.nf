
process MEAN_TCR_PMHC_PAE {
  label "process_local"
  conda "envs/env.yaml"

  publishDir "${params.outdir}/triad/features", mode: 'copy'

  input:
  path triad_pq

  output:
  path("*.parquet")

  script:
  """
  #!/usr/bin/env python
  from tcrtrifold.feat_extract import extract_mean_tcr_pmhc_pae
  from mdaf3.FeatureExtraction import split_apply_combine
  import polars as pl
  from pathlib import Path

  df = pl.read_parquet("${triad_pq}")
  df = split_apply_combine(
    df,
    extract_mean_tcr_pmhc_pae,
    Path("${params.outdir}/triad/inference"),
    chunksize=15,
  )
  df = df.select("job_name", "mean_p_tcr_pae")
  df.write_parquet("${triad_pq.getBaseName()}.mean_tcr_pmhc_pae.parquet")
  """
}


