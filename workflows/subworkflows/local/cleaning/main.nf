


// process CLEAN_CRESTA_WINDOWED {
//   label "tcrtrifold_local"

//   publishDir(
//     path: {"${params.outdir}/triad/staged"},
//     pattern: "*triad*",
//     mode: 'copy'
//   )
//   publishDir(
//     path: {"${params.outdir}/pmhc/staged"},
//     pattern: "*pmhc*",
//     mode: 'copy'
//   )


//   input:
//   path cresta
//   path peptide_correction_csv

//   output:
//   path("*triad*.parquet"), emit: triad
//   path("*pmhc*.parquet"), emit: pmhc

//   script:
//   """
//   clean_cresta_w.py \\
//     --raw_csv_path ${cresta} \\
//     --peptide_correction_csv ${peptide_correction_csv} \\
//     -ot cresta_triad.cleaned.parquet \\
//     -op cresta_pmhc.cleaned.parquet 
//   """
// }

process CLEAN_CRESTA {
  label "tcrtrifold_local"

  // publishDir(
  //   path: {"${params.outdir}/triad/staged"},
  //   pattern: "*triad*",
  //   mode: 'copy'
  // )
  // publishDir(
  //   path: {"${params.outdir}/pmhc/staged"},
  //   pattern: "*pmhc*",
  //   mode: 'copy'
  // )


  input:
  path cresta

  output:
  path("*triad*.parquet"), emit: triad_parquet
  path("*pmhc*.parquet"), emit: pmhc_parquet

  script:
  """
  clean_cresta.py \\
    --raw_csv_path ${cresta} \\
    -ot cresta_triad.cleaned.parquet \\
    -op cresta_pmhc.cleaned.parquet 
  """
}

process CLEAN_IEDB_I {
  label "tcrtrifold_heavy"

  input:
  path iedb_I

  output:
  path("*triad*.parquet"), emit: triad_parquet
  path("*pmhc*.parquet"), emit: pmhc_parquet

  script:
  """
  clean_iedb_I.py \\
    --raw_csv_path ${iedb_I} \\
    -ot iedb_I_triad.cleaned.parquet \\
    -op iedb_I_pmhc.cleaned.parquet 
  """
}


process CLEAN_IEDB_II {
  label "tcrtrifold_heavy"

  input:
  path iedb_II

  output:
  path("*triad*.parquet"), emit: triad_parquet
  path("*pmhc*.parquet"), emit: pmhc_parquet

  script:
  """
  clean_iedb_II.py \\
    --raw_csv_path ${iedb_II} \\
    -ot iedb_II_triad.cleaned.parquet \\
    -op iedb_II_pmhc.cleaned.parquet 
  """
}



process CLEAN_PDB {
  label "tcrtrifold_local"
  // don't publish anything, they might be filtered out in next step

  input:
  path pdb_rep
  path pdb_stcr

  output:
  path("*triad*.parquet"), emit: triad_parquet
  path("*pmhc*.parquet"), emit: pmhc_parquet

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

  input:
  path pdb_pq

  output:
  path("*.pdb")

  script:
  """





  
  format_true_pdbs.py \\
    --pdb_parquet ${pdb_pq} \\
    --output_dir .
  """
}

process VALIDATION_STANDARDIZE {
    label "tcrtrifold_local"

  input:
    path pdb_pq
  output:
    path("*triad*.parquet"), emit: triad
    path("*pmhc*.parquet"), emit: pmhc

  script:
  """
  validation_standardize.py \\
    --pdb_parquet ${pdb_pq} \\
    -ot pdb_validation_triad.cleaned.parquet \\
    -op pdb_validation_pmhc.cleaned.parquet 
  """
}


// workflow MSA_WORKFLOW {
//     take:
//     meta_fasta

//     main:
//     FILT_FORMAT_MSA(meta_fasta)
//     RUN_MSA(FILT_FORMAT_MSA.out)
//     STORE_MSA(RUN_MSA.out)

//     emit:
//     new_msa = STORE_MSA
// }


// workflow CLEAN_PDB_VALIDATION_WORKFLOW {

//     main:

//     CLEAN_PDB(Channel.fromPath("${params.input}/table_S1_structure_benchmark_complexes.csv"),
//         Channel.fromPath("${params.input}/raw/db_summary.dat"))

//     VALIDATION_STANDARDIZE(CLEAN_PDB.out.triad)

//     FORMAT_TRUE_PDBS(VALIDATION_STANDARDIZE.out.triad)

//     emit:
//     triad = FORMAT_TRUE_PDBS.out
//     pmhc = VALIDATION_STANDARDIZE.out.pmhc

// }

// workflow CLEAN_CRESTA_WINDOWED_WORKFLOW {
//     take:
//     raw_data_path

//     main:
//     CLEAN_CRESTA_WINDOWED(Channel.fromPath("${params.input}/cresta.csv"), Channel.fromPath("${params.input}/peptide_corrections.csv"))
//     emit:
//     triad = CLEAN_CRESTA_WINDOWED.out.triad
//     pmhc = CLEAN_CRESTA_WINDOWED.out.pmhc
// }

// workflow {

//   CLEAN_IEDB_II(Channel.fromPath("${params.data_dir}/${params.dset_name}/raw/immrep_IEDB.csv"))

// }


// workflow {
//   CLEAN_IEDB_I(Channel.fromPath("${params.data_dir}/${params.dset_name}/raw/immrep_IEDB.csv"))
// }

// workflow {

//   CLEAN_CRESTA(Channel.fromPath("${params.data_dir}/${params.dset_name}/raw/cresta.csv"))
// }

