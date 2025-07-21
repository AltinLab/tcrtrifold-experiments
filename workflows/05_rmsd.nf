process COMPUTE_RMSD {
  label "tcrtrifold_local"

  publishDir "${params.data_dir}/pdb/triad/staged", mode: 'copy'

  input:
  path pdb_pq
  path cleaned_pdb_dir
  path inference_dir
  path true_tcrdock_pq
  path pred_tcrdock_pq

  output:
  path("*triad*.parquet"), emit: triad

  script:
  """
  rmsd.py \\
    --input_parquet ${pdb_pq} \\
    --cleaned_pdbs ${cleaned_pdb_dir} \\
    --inference_type ${params.inf_type} \\
    --inference_dir ${inference_dir} \\
    --true_tcrdock_pq ${true_tcrdock_pq} \\
    --pred_tcrdock_pq ${pred_tcrdock_pq} \\
    -o ${pdb_pq.getSimpleName()}.${params.inf_type}_rmsd.parquet
  """
}


workflow {

    pq = Channel.fromPath("data/pdb/triad/staged/*.cleaned*.parquet")
    cleaned_pdbs = Channel.fromPath("data/pdb/triad/cleaned_pdb")
    true_tcrdock_pq = Channel.fromPath("data/pdb/triad/staged/*.true_tcrdock*.parquet")
    pred_tcrdock_pq = Channel.fromPath("data/pdb/triad/staged/*.${params.inf_type}_tcrdock*.parquet")

    if (params.inf_type == "af3") {
        inf_dir = Channel.fromPath("data/pdb/triad/inference")
    }
    else {
        inf_dir = Channel.fromPath("data/pdb/triad/predictions")
    }
    

    COMPUTE_RMSD(pq, cleaned_pdbs, inf_dir, true_tcrdock_pq, pred_tcrdock_pq)

}