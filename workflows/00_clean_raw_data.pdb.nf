/*
Hardcoded pipelines for cleaning raw data files from each dataset, getting them into a uniform format.
*/

process CLEAN_PDB {
  label "tcrtrifold_local"

  publishDir(
    path: {"${params.data_dir}/pdb/triad/staged"},
    pattern: "*triad*",
    mode: 'copy'
  )
  publishDir(
    path: {"${params.data_dir}/pdb/pmhc/staged"},
    pattern: "*pmhc*",
    mode: 'copy'
  )


  input:
  path pdb_rep
  path pdb_stcr

  output:
  path("*triad*.parquet"), emit: triad
  path("*pmhc*.parquet"), emit: pmhc

  script:
  """
  clean_pdb.py \
    --raw_csv_path ${pdb_rep} \\
    --raw_stcr_path ${pdb_stcr} \\
    --imgt_hla_path ${params.imgt_hla_path} \\
    -ot pdb_triad.cleaned.parquet \\
    -op pdb_pmhc.cleaned.parquet 
  """
}

process FORMAT_TRUE_PDBS {
  label "tcrtrifold_local"

  publishDir "${params.data_dir}/pdb/triad/cleaned_pdb", mode: 'copy'

  input:
  path pdb_pq

  output:
  path("*.pdb")

  script:
  """
  format_true_pdbs.py \
    --pdb_parquet ${pdb_pq} \\
    --output_dir .
  """
}


workflow {

  CLEAN_PDB(Channel.fromPath("data/pdb/raw/table_S1_structure_benchmark_complexes.csv"),
    Channel.fromPath("data/pdb/raw/db_summary.dat"))

  FORMAT_TRUE_PDBS(CLEAN_PDB.out.triad)

}