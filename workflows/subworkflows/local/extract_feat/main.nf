
process EXTRACT_TRIAD_CONF_FEAT {
  label "tcrtrifold_local"

  input:
  path triad_pq
  val inference_dir
  val inf_type

  output:
  path("*.parquet")

  script:
  """




  extract_conf_feat.py \\
      --input_parquet ${triad_pq} \\
      --inference_dir ${inference_dir} \\
      --inference_type ${inf_type} \\
      --output_path "${triad_pq.getSimpleName()}.conf_${inf_type}.parquet"
  """
}


process EXTRACT_PMHC_CONF_FEAT {
  label "tcrtrifold_local"

  input:
  path pmhc_pq
  val inference_dir
  val inf_type

  output:
  path("*.parquet")

  script:
  """
  extract_pmhc_conf_feat.py \\
      --input_parquet ${pmhc_pq} \\
      --inference_dir ${inference_dir} \\
      --inference_type ${inf_type} \\
      --output_path "${pmhc_pq.getSimpleName()}.conf_${inf_type}.parquet"
  """
}

process EXTRACT_TRIAD_SUMMARY_METRICS {
  label "tcrtrifold_local"


  input:
  path triad_pq
  val inference_dir
  val inf_type

  output:
  path("*.parquet")

  script:
  """
  extract_conf_feat.py \\
      --input_parquet ${triad_pq} \\
      --inference_dir ${inference_dir} \\
      --inference_type ${inf_type} \\
      --summary_only \\
      --output_path "${triad_pq.getSimpleName()}.conf_${inf_type}.parquet"
  """
}

process TCRDOCK_GEOM_FROM_PDB {
    label "tcrdock"

    input:
    path input_pq
    path topology_files

    output:
    path("*.parquet")

    script:
    """

    tcrdock_geom.py \\
        --input_parquet ${input_pq} \\
        --topology_path . \\
        --from_true_struct \\
        -o "${input_pq.getSimpleName()}.true_tcrdock.parquet"
    """
    }     
    


process TCRDOCK_GEOM_FROM_INFERENCE {
    label "tcrdock"

    input:
    path input_pq
    val inf_dir
    val inf_type

    output:
    path("*.parquet")

    script:
    """
    
    tcrdock_geom.py \\
        --input_parquet ${input_pq} \\
        --topology_path ${inf_dir} \\
        --inference_type ${inf_type} \\
        -o "${input_pq.getSimpleName()}.${inf_type}_tcrdock.parquet"
    """
}

process COMPUTE_RMSD {
  label "tcrtrifold_local"

  input:
  path triad_parquet
  path true_tcrdock_pq
  path pred_tcrdock_pq
  path cleaned_pdbfiles
  path inference_dir
  val inf_type

  output:
  path("*triad*.parquet"), emit: triad

  script:
  """

  rmsd.py \\
    --input_parquet ${triad_parquet} \\
    --cleaned_pdbs . \\
    --inference_type ${inf_type} \\
    --inference_dir ${inference_dir} \\
    --true_tcrdock_pq ${true_tcrdock_pq} \\
    --pred_tcrdock_pq ${pred_tcrdock_pq} \\
    -o ${true_tcrdock_pq.getSimpleName()}.${inf_type}_rmsd.parquet
  """
}



process PARQUET_TO_FASTA {
  label "tcrtrifold_local"
  
  input:
      tuple val(meta), path(parquet)
  
  output:
      tuple val(meta), path("*.fasta")
  
  script:
  """
  #!/usr/bin/env python

  import polars as pl
  
  df = pl.read_parquet("${parquet}")
  
  with open("${parquet.getSimpleName()}.fasta", "w") as f:
      for row in df.iter_rows(named=True):
          f.write(f">{row['seq']}\\n{row['seq']}\\n")
  """
}


process CLUSTER_FASTA {
    label "process_local"
    conda "envs/mmseqs2.yaml"

    input:
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path("*.fasta")

    script:
    """
    mmseqs createdb ${fasta} DB && \\
    mmseqs cluster DB DB_clu tmp \\
        --min-seq-id 0.8 \\
        --cov-mode 1 && \\
    mmseqs createsubdb DB_clu DB DB_clu_rep && \\
    mmseqs convert2fasta DB_clu_rep "${fasta.getSimpleName()}.clust.fasta"
    """
}



process FASTA_TO_NEFF {
  label "tcrtrifold_local"
  
  input:
      tuple val(meta), path(fasta)
  
  output:
      path("*.parquet")
  
  script:
  """
  #!/usr/bin/env python

  import polars as pl
  from tcrtrifold.utils import fasta_to_polars
  
  df = fasta_to_polars("${fasta}")
  
  out_df = pl.DataFrame(
      {
          "seq": "${meta.id}",
          "protein_chain": "${meta.protein_chain}",
          "neff": df.height,
      }
  )

  out_df.write_parquet("${fasta.getSimpleName()}.parquet")
  """
}


workflow NEFF_WORKFLOW {
    take:
    chain_msa_fasta

    main:
    CLUSTER_FASTA(chain_msa_fasta)
    out_ch = FASTA_TO_NEFF(CLUSTER_FASTA.out)

    emit:
    out_ch
}