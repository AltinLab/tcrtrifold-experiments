
// process EXCLUDE_VALIDATION {
//   label "tcrtrifold_local"

//   input:
//   path triad_pq
//   path validation_exclusions

//   output:
//   path("*triad*.parquet")

//   script:
//   """
//   exclude_validation.py \\
//     --triad_parquet ${triad_pq} \\
//     --validation_exclusions ${validation_exclusions} \\
//     -ot ${triad_pq.getSimpleName()}.cleaned.parquet
//   """
// }

process SPLIT_TRIAD_INTO_CHAINS {
  label "tcrtrifold_local"
  
  input:
  path triad_parquet

  output:
  path "*hashed_seq*.parquet", emit: triad_pq
  path "*peptide*.parquet", emit: peptide_pq
  path "*mhc_1*.parquet", emit: mhc_1_pq
  path "*mhc_2*.parquet", emit: mhc_2_pq
  path "*tcr_1*.parquet", emit: tcr_1_pq
  path "*tcr_2*.parquet", emit: tcr_2_pq

  script:
  """

  split_triad_into_chains.py \\
    --triad_parquet ${triad_parquet} \\
    --output_triad ${triad_parquet.getSimpleName()}.hashed_seq.parquet \\
    --output_peptide ${triad_parquet.getSimpleName()}_peptide.parquet \\
    --output_mhc_1 ${triad_parquet.getSimpleName()}_mhc_1.parquet \\
    --output_mhc_2 ${triad_parquet.getSimpleName()}_mhc_2.parquet \\
    --output_tcr_1 ${triad_parquet.getSimpleName()}_tcr_1.parquet \\
    --output_tcr_2 ${triad_parquet.getSimpleName()}_tcr_2.parquet
  
  """

}

process PARQUET_TO_FASTA {
  label "tcrtrifold_local"

  input:
      path parquet
  
  output:
      path "*.fasta", emit: fasta
  
  script:
  """
  #!/usr/bin/env python

  import polars as pl
  
  df = pl.read_parquet("${parquet}").filter(pl.col("seq").is_not_null())
  
  with open("${parquet.getSimpleName()}.fasta", "w") as f:
      for row in df.iter_rows(named=True):
          f.write(f">{row['name']}\\n{row['seq']}\\n")
  """
}

process MMSEQS_QUERY_TARGET_IDENT {
  label "mmseqs_heavy"

  input:
  path query_fasta
  path target_fasta

  output:
  path "*.tsv", emit: ident_tsv

  script:
  """
  mmseqs createdb ${query_fasta} query_db
  mmseqs createdb ${target_fasta} target_db
  mmseqs search query_db target_db ident_db tmpdir --alignment-mode 3 -s 7.5 --min-seq-id 0.6
  mmseqs convertalis query_db target_db ident_db "${query_fasta.getSimpleName()}_${target_fasta.getSimpleName()}.tsv" \\
    --format-output "query,target,pident,alnlen,qcov,tcov,evalue,bits,qstart,qend,tstart,tend"
  """
}

process MMSEQS_QUERY_TARGET_IDENT_SHORT {
  label "mmseqs_heavy"

  input:
  path query_fasta
  path target_fasta

  output:
  path "*.tsv", emit: ident_tsv

  script:
  """
  mmseqs createdb ${query_fasta} query_db
  mmseqs createdb ${target_fasta} target_db
  # see: https://github.com/soedinglab/MMseqs2/issues/125 for why spaced-kmer-mode arg is needed
  mmseqs search query_db target_db ident_db tmpdir --alignment-mode 3 -s 7.5 --min-seq-id 0.6 --spaced-kmer-mode 0
  mmseqs convertalis query_db target_db ident_db "${query_fasta.getSimpleName()}_${target_fasta.getSimpleName()}.tsv" \\
    --format-output "query,target,pident,alnlen,qcov,tcov,evalue,bits,qstart,qend,tstart,tend"
  """
}

process FILTER_ANNOTATE_TRIAD_FROM_IDENT {
  label "tcrtrifold_local"

  input:
  path query_triad_pq
  path curr_pmhc_parquet
  path target_triad_pq
  path peptide_ident_tsv
  path mhc_1_ident_tsv
  path mhc_2_ident_tsv
  path tcr_1_ident_tsv
  path tcr_2_ident_tsv

  output:
  path "*annotated*.parquet", emit: annot_triad_parquet
  path "*excluded_triad*.parquet", emit: excluded_triad_parquet
  path "*remaining*.parquet", emit: remaining_pmhc_parquet

  script:
  def base_name = query_triad_pq.getSimpleName().substring(0, query_triad_pq.getSimpleName().indexOf('_triad'))
  """
  filter_annotate_triad_from_ident.py \\
    --query_triad_parquet ${query_triad_pq} \\
    --pmhc_parquet ${curr_pmhc_parquet} \\
    --validation_triad_parquet ${target_triad_pq} \\
    --peptide_ident_tsv ${peptide_ident_tsv} \\
    --mhc_1_ident_tsv ${mhc_1_ident_tsv} \\
    --mhc_2_ident_tsv ${mhc_2_ident_tsv} \\
    --tcr_1_ident_tsv ${tcr_1_ident_tsv} \\
    --tcr_2_ident_tsv ${tcr_2_ident_tsv} \\
    --output_annotated_triad_parquet ${query_triad_pq.getSimpleName()}.annotated.parquet \\
    --output_excluded_triad_parquet ${query_triad_pq.getSimpleName()}.excluded_triad.parquet \\
    --output_remaining_pmhc_parquet ${base_name}_pmhc.remaining.parquet
  """
}

process FILTER_ANNOTATE_TRIAD_FROM_PDB {
  label "tcrtrifold_local"

  input:
  path query_triad_pq
  path pmhc_parquet
  path blast_mhc_1_tsv
  path blast_mhc_2_tsv
  path blast_tcr_1_tsv
  path blast_tcr_2_tsv

  output:
  path "*annotated*.parquet", emit: annot_triad_parquet
  path "*excluded_triad*.parquet", emit: excluded_triad_parquet
  path "*remaining*.parquet", emit: remaining_pmhc_parquet

  script:
  def base_name = query_triad_pq.getSimpleName().substring(0, query_triad_pq.getSimpleName().indexOf('_triad'))
  """


  filter_annotate_triad_from_pdb.py \\
    --query_triad_parquet ${query_triad_pq} \\
    --pmhc_parquet ${pmhc_parquet} \\
    --blast_mhc_1_tsv ${blast_mhc_1_tsv} \\
    --blast_mhc_2_tsv ${blast_mhc_2_tsv} \\
    --blast_tcr_1_tsv ${blast_tcr_1_tsv} \\
    --blast_tcr_2_tsv ${blast_tcr_2_tsv} \\
    --output_annotated_triad_parquet ${query_triad_pq.getSimpleName()}.annotated.parquet \\
    --output_excluded_triad_parquet ${query_triad_pq.getSimpleName()}.excluded_triad.parquet \\
    --output_remaining_pmhc_parquet ${base_name}_pmhc.remaining.parquet \\
    --exclude_af3_training_only
  """
}

process ANNOTATE_TRIAD_FROM_PDB {
  label "tcrtrifold_local"

  input:
  path query_triad_pq
  path blast_mhc_1_tsv
  path blast_mhc_2_tsv
  path blast_tcr_1_tsv
  path blast_tcr_2_tsv

  output:
  path "*annotated*.parquet", emit: annot_triad_parquet

  script:
  """


  filter_annotate_triad_from_pdb.py \\
    --query_triad_parquet ${query_triad_pq} \\
    --blast_mhc_1_tsv ${blast_mhc_1_tsv} \\
    --blast_mhc_2_tsv ${blast_mhc_2_tsv} \\
    --blast_tcr_1_tsv ${blast_tcr_1_tsv} \\
    --blast_tcr_2_tsv ${blast_tcr_2_tsv} \\
    --output_annotated_triad_parquet ${query_triad_pq.getSimpleName()}.annotated.parquet \\
    --exclude_af3_training_only \\
    --annotate_only
  """
}

process BLAST_PARQUET_WEBSERVER {
  label "tcrtrifold_local"

  input:
  path seq_pq

  output:
  path "*.tsv", emit: blast_tsv

  script:
  """
  blast_parquet_webserver.py \\
    --seq_parquet ${seq_pq} \\
    --output_blast_tsv ${seq_pq.getSimpleName()}.blast.tsv
  """
}