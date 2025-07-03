/*
Hardcoded pipelines for cleaning raw data files from each dataset, getting them into a uniform format.
*/

process CLEAN_CRESTA {
  label "tcrtrifold_local"

  publishDir(
    path: {"${params.data_dir}/cresta/triad/staged"},
    pattern: "*triad*",
    mode: 'copy'
  )
  publishDir(
    path: {"${params.data_dir}/cresta/pmhc/staged"},
    pattern: "*pmhc*",
    mode: 'copy'
  )


  input:
  path cresta

  output:
  path("*triad*.parquet")
  path("*pmhc*.parquet")

  script:
  """
  clean_cresta.py \\
    --raw_csv_path ${cresta} \\
    -ot cresta_triad.cleaned.parquet \\
    -op cresta_pmhc.cleaned.parquet 
  """
}


workflow {

  CLEAN_CRESTA(Channel.fromPath("data/cresta/raw/cresta.csv"))
}