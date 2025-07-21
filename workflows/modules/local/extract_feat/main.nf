
process EXTRACT_TRIAD_CONF_FEAT {
  label "process_local"
  conda "envs/env.yaml"

  publishDir "${params.outdir}/triad/staged", mode: 'copy'

  input:
  path triad_pq
  path inference_dir

  output:
  path("*.parquet")

  script:
  """
  extract_conf_feat.py \\
      --input_parquet ${triad_pq} \\
      --inference_dir ${inference_dir} \\
      --inference_type ${params.inf_type} \\
      --output_path "${triad_pq.getSimpleName()}.conf.parquet"
  """
}


