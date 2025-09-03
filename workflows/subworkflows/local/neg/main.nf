process GEN_NEGATIVES {
  label "tcrtrifold_very_heavy"


  input:
  path base_df
  val supp_dfs
  output:
  path("*neg*.parquet"), emit: triad_negatives
  path("*discard*.parquet"), optional: true, emit: discard_df

  script:
  """
  gen_negatives.py \\
    --base_df ${base_df} \\
    --supp_dfs ${supp_dfs.join(",")} \\
    --discard_path "${base_df.getSimpleName()}.discard.parquet" \\
    --output_path "${base_df.getSimpleName()}.neg.parquet"
  """
}

process PERFORM_WINDOW {
  label "tcrtrifold_local"

  input:
  path neg_df

  output:
  path("*triad*.parquet"), emit: triad_df
  path("*pmhc*.parquet"), emit: pmhc_df

  script:
  def base_name = neg_df.getSimpleName().substring(0, neg_df.getSimpleName().indexOf('_triad'))
  """
  window_peptides.py \\
    --neg_df ${neg_df} \\
    --output_pmhc_path "${base_name}_pmhc_w.cleaned.parquet" \\
    --output_triad_path "${neg_df.getSimpleName()}_w.neg.parquet"
  """
}
