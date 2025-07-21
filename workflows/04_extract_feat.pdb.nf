process COMPUTE_RMSD {
  label "tcrtrifold_local"

  publishDir "${params.data_dir}/pdb/triad/staged", mode: 'copy'

  input:
  path pdb_pq
  path cleaned_pdb_dir
  path inference_dir

  output:
  path("*triad*.parquet"), emit: triad

  script:
  """
  compute_rmsd.py \
    --input_parquet ${pdb_pq} \\
    --cleaned_pdbs ${cleaned_pdb_dir} \\
    --inference_dir ${inference_dir} \\
    -o ${pdb_pq.getSimpleName()}.rmsd.parquet
  """
}


workflow {

    clean_pq = Channel.fromPath("data/pdb/triad/staged/*.cleaned*.parquet")

    neg_pq = Channel.fromPath("data/pdb/triad/staged/*.neg*.parquet")
    inf_dir = Channel.fromPath("data/pdb/triad/inference")
    cleaned_pdbs = Channel.fromPath("data/pdb/triad/cleaned_pdb")

    COMPUTE_RMSD(clean_pq, cleaned_pdbs, inf_dir)

}