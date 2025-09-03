/*
Parameters:
- input_replication: Path to the replication CSV file
- input_stcrdb_summary: Path to the STCRDB summary CSV file
*/

// in current version, new output syntax is in preview
nextflow.preview.output = true


include { MSA_WORKFLOW } from './subworkflows/tgen/af3'
include { SEQ_LIST_TO_FASTA } from './modules/tgen/af3'
include { CLEAN_PDB; FORMAT_TRUE_PDBS; VALIDATION_STANDARDIZE } from './subworkflows/local/cleaning'
include { EXCLUDE_AF3_TRAINING_DATA } from './subworkflows/local/clustering'
include { GEN_NEGATIVES; PERFORM_WINDOW } from './subworkflows/local/neg'
include { MSA_FROM_TRIAD_PARQUET; 
        UNBATCHED_INFERENCE_FROM_TRIAD_PARQUET as UNBATCHED_INFERENCE_ALL_PDB; 
        UNBATCHED_INFERENCE_FROM_TRIAD_PARQUET as UNBATCHED_INFERENCE_VALIDATION;
        UNBATCHED_INFERENCE_FROM_PMHC_PARQUET;
        NOOP_DEP as DEPEND_TRIAD_ALL_ON_INFERENCE; 
        NOOP_DEP as DEPEND_TRIAD_ALL_ON_BOLTZ; 
        NOOP_DEP as DEPEND_TRIAD_VALIDATION_ON_INFERENCE;
        NOOP_DEP as DEPEND_VALIDATION_ON_ALL_INFERENCE;
        NOOP_DEP as DEPEND_PMHC_ON_INFERENCE  } from './subworkflows/local/af3_adapter'
include { BOLTZ_FROM_TRIAD_PARQUET } from './subworkflows/local/boltz'
include { EXTRACT_TRIAD_SUMMARY_METRICS;
        EXTRACT_TRIAD_SUMMARY_METRICS as EXTRACT_TRIAD_SUMMARY_METRICS_BOLTZ;
        EXTRACT_TRIAD_CONF_FEAT;
        TCRDOCK_GEOM_FROM_PDB; 
        TCRDOCK_GEOM_FROM_INFERENCE as TCRDOCK_GEOM_FROM_AF3_INFERENCE;
        TCRDOCK_GEOM_FROM_INFERENCE as TCRDOCK_GEOM_FROM_BOLTZ_INFERENCE;
        TCRDOCK_GEOM_FROM_INFERENCE as TCRDOCK_GEOM_VALIDATION;
        COMPUTE_RMSD as COMPUTE_RMSD_AF3;
        COMPUTE_RMSD as COMPUTE_RMSD_BOLTZ } from './subworkflows/local/extract_feat'

workflow {

    main:

    rep_csv = Channel.fromPath(params.input_replication)
    stcr_csv = Channel.fromPath(params.input_stcr)

    triad_inf_dir = file("$workflow.outputDir/triad/inference").toUriString()
    pmhc_inf_dir = file("$workflow.outputDir/pmhc/inference").toUriString()
    triad_boltz_dir = file("$workflow.outputDir/triad/predictions").toUriString()
    pmhc_boltz_dir = file("$workflow.outputDir/pmhc/predictions").toUriString()

    CLEAN_PDB(rep_csv, stcr_csv)
    triad_cleaned = CLEAN_PDB.out.triad_parquet
    pmhc_cleaned = CLEAN_PDB.out.pmhc_parquet

    cleaned_pdbfiles = FORMAT_TRUE_PDBS(triad_cleaned)

    // ALL PDB inference

    MSA_FROM_TRIAD_PARQUET(triad_cleaned)
    triad_msa_done_token = MSA_FROM_TRIAD_PARQUET.out.new_msa_list

    UNBATCHED_INFERENCE_ALL_PDB(triad_cleaned, triad_inf_dir, triad_msa_done_token)
    triad_all_meta_inf = UNBATCHED_INFERENCE_ALL_PDB.out.new_meta_inf
    triad_cleaned_af3 = DEPEND_TRIAD_ALL_ON_INFERENCE(triad_cleaned, triad_all_meta_inf.toList())

    BOLTZ_FROM_TRIAD_PARQUET(triad_cleaned, triad_boltz_dir)
    triad_all_meta_boltz = BOLTZ_FROM_TRIAD_PARQUET.out.new_meta_inf
    triad_cleaned_boltz = DEPEND_TRIAD_ALL_ON_BOLTZ(triad_cleaned, triad_all_meta_boltz.toList())

    triad_tcrdock_pdb = TCRDOCK_GEOM_FROM_PDB(triad_cleaned, cleaned_pdbfiles)
    triad_tcrdock_af3 = TCRDOCK_GEOM_FROM_AF3_INFERENCE(triad_cleaned_af3, triad_inf_dir, Channel.value("af3"))
    triad_tcrdock_boltz = TCRDOCK_GEOM_FROM_BOLTZ_INFERENCE(triad_cleaned_boltz, triad_boltz_dir, Channel.value("boltz"))

    COMPUTE_RMSD_AF3(triad_cleaned, triad_tcrdock_pdb, triad_tcrdock_af3, cleaned_pdbfiles, triad_inf_dir, Channel.value("af3"))
    triad_rmsd_af3 = COMPUTE_RMSD_AF3.out.triad

    COMPUTE_RMSD_BOLTZ(triad_cleaned, triad_tcrdock_pdb, triad_tcrdock_boltz, cleaned_pdbfiles, triad_boltz_dir, Channel.value("boltz"))
    triad_rmsd_boltz = COMPUTE_RMSD_BOLTZ.out.triad

    triad_all_conf = EXTRACT_TRIAD_SUMMARY_METRICS(triad_cleaned_af3, triad_inf_dir, Channel.value("af3"))

    triad_all_conf_boltz = EXTRACT_TRIAD_SUMMARY_METRICS_BOLTZ(triad_cleaned_boltz, triad_boltz_dir, Channel.value("boltz"))

    // VALIDATION SET

    VALIDATION_STANDARDIZE(triad_cleaned)
    triad_validation_cleaned = VALIDATION_STANDARDIZE.out.triad
    pmhc_validation_cleaned = VALIDATION_STANDARDIZE.out.pmhc


    EXCLUDE_AF3_TRAINING_DATA(triad_validation_cleaned, pmhc_validation_cleaned)
    triad_validation_cleaned = EXCLUDE_AF3_TRAINING_DATA.out.annot_triad_parquet
    triad_validation_excluded_triad = EXCLUDE_AF3_TRAINING_DATA.out.excluded_triad_parquet
    pmhc_validation_cleaned = EXCLUDE_AF3_TRAINING_DATA.out.remaining_pmhc_parquet

    GEN_NEGATIVES(triad_validation_cleaned, triad_validation_cleaned.toList())
    triad_validation_negatives = GEN_NEGATIVES.out.triad_negatives

    UNBATCHED_INFERENCE_VALIDATION(triad_validation_negatives, triad_inf_dir, triad_all_meta_inf.toList())
    triad_validation_meta_inf = UNBATCHED_INFERENCE_VALIDATION.out.new_meta_inf
    triad_validation_negatives = DEPEND_TRIAD_VALIDATION_ON_INFERENCE(triad_validation_negatives, triad_validation_meta_inf.toList())

    triad_conf_validation = EXTRACT_TRIAD_CONF_FEAT(triad_validation_negatives, triad_inf_dir, Channel.value("af3"))

    triad_tcrdock_validation = TCRDOCK_GEOM_VALIDATION(triad_validation_negatives, triad_inf_dir, Channel.value("af3"))

    UNBATCHED_INFERENCE_FROM_PMHC_PARQUET(pmhc_cleaned, pmhc_inf_dir, triad_msa_done_token)
    pmhc_validation_meta_inf = UNBATCHED_INFERENCE_FROM_PMHC_PARQUET.out.new_meta_inf
    pmhc_cleaned = DEPEND_PMHC_ON_INFERENCE(pmhc_cleaned, pmhc_validation_meta_inf.toList())

    publish:
    triad_cleaned = triad_cleaned
    pmhc_cleaned = pmhc_cleaned
    cleaned_pdbfiles = cleaned_pdbfiles
    incomplete_negative_log = GEN_NEGATIVES.out.discard_df

    triad_all_meta_inf = triad_all_meta_inf
    triad_all_meta_boltz = triad_all_meta_boltz

    triad_tcrdock_pdb = triad_tcrdock_pdb
    triad_tcrdock_af3 = triad_tcrdock_af3
    triad_tcrdock_boltz = triad_tcrdock_boltz

    triad_rmsd_af3 = triad_rmsd_af3
    triad_rmsd_boltz = triad_rmsd_boltz

    triad_all_conf = triad_all_conf
    triad_all_conf_boltz = triad_all_conf_boltz

    triad_validation_cleaned = triad_validation_cleaned
    triad_validation_excluded_triad = triad_validation_excluded_triad
    triad_validation_negatives = triad_validation_negatives
    triad_validation_meta_inf = triad_validation_meta_inf
    pmhc_validation_meta_inf = pmhc_validation_meta_inf

    triad_conf_validation = triad_conf_validation
    triad_tcrdock_validation = triad_tcrdock_validation

}

output {
    triad_cleaned {
        path "triad/staged"
    }
    pmhc_cleaned {
        path "pmhc/staged"
    }
    cleaned_pdbfiles {
        path "triad/cleaned_pdb"
    }
    incomplete_negative_log {
        path "triad/staged"
    }

    triad_all_meta_inf {
        path "triad/inference"
    }
    triad_all_meta_boltz {
        path "triad/predictions"
    }
    
    triad_tcrdock_pdb {
        path "triad/staged"
    }
    triad_tcrdock_af3 {
        path "triad/staged"
    }
    triad_tcrdock_boltz {
        path "triad/staged"
    }

    triad_rmsd_af3 {
        path "triad/staged"
    }
    triad_rmsd_boltz {
        path "triad/staged"
    }

    triad_all_conf {
        path "triad/staged"
    }
    triad_all_conf_boltz {
        path "triad/staged"
    }
    
    triad_validation_cleaned {
        path "triad/staged"
    }
    triad_validation_excluded_triad {
        path "triad/staged"
    }
    triad_validation_negatives {
        path "triad/staged"
    }
    triad_validation_meta_inf {
        path "triad/inference"
    }
    pmhc_validation_meta_inf {
        path "pmhc/inference"
    }
    triad_conf_validation {
        path "triad/staged"
    }
    triad_tcrdock_validation {
        path "triad/staged"
    }

}
