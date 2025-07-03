/*
Hardcoded pipelines for cleaning raw data files from each dataset, getting them into a uniform format.
*/

process CLEAN_IEDB_II {
  label "tcrtrifold_heavy"

  publishDir(
      path: {"${params.data_dir}/${params.dset_name}/triad/staged"},
      pattern: "*triad*",
      mode: 'copy'
  )
  publishDir(
      path: {"${params.data_dir}/${params.dset_name}/pmhc/staged"},
      pattern: "*pmhc*",
      mode: 'copy'
  )

  input:
  path iedb_II

  output:
  path("*.parquet")

  script:
  """
  clean_iedb_II.py \\
    --raw_csv_path ${iedb_II} \\
    -ot ${params.dset_name}_triad.cleaned.parquet \\
    -op ${params.dset_name}_pmhc.cleaned.parquet 
  """
}


workflow {

  CLEAN_IEDB_II(Channel.fromPath("${params.data_dir}/${params.dset_name}/raw/immrep_IEDB.csv"))

}