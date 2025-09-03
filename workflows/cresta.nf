/*
Parameters:
- input: Path to the CRESTA raw data CSV file
*/

// in current version, new output syntax is in preview
nextflow.preview.output = true


include { MSA_WORKFLOW } from './subworkflows/tgen/af3'
include { SEQ_LIST_TO_FASTA } from './modules/tgen/af3'
include { CLEAN_CRESTA } from './subworkflows/local/cleaning'
include { EXCLUDE_AF3_TRAINING_DATA } from './subworkflows/local/clustering'
include { GEN_NEGATIVES; PERFORM_WINDOW } from './subworkflows/local/neg'
include { MSA_FROM_TRIAD_PARQUET; 
        UNBATCHED_INFERENCE_FROM_TRIAD_PARQUET as UNBATCHED_INFERENCE_FROM_UNWIN_TRIAD_PARQUET; 
        UNBATCHED_INFERENCE_FROM_PMHC_PARQUET as UNBATCHED_INFERENCE_FROM_UNWIN_PMHC_PARQUET;
        NOOP_DEP as DEPEND_TRIAD_ON_UNWIN_INFERENCE; 
        NOOP_DEP as DEPEND_PMHC_ON_UNWIN_INFERENCE;
        UNBATCHED_INFERENCE_FROM_TRIAD_PARQUET; 
        UNBATCHED_INFERENCE_FROM_PMHC_PARQUET;
        NOOP_DEP as DEPEND_TRIAD_ON_INFERENCE; 
        NOOP_DEP as DEPEND_PMHC_ON_INFERENCE  } from './subworkflows/local/af3_adapter'
include { EXTRACT_TRIAD_CONF_FEAT; 
    EXTRACT_TRIAD_CONF_FEAT as EXTRACT_UNWIN_TRIAD_CONF_FEAT; 
    EXTRACT_PMHC_CONF_FEAT; 
    TCRDOCK_GEOM_FROM_INFERENCE } from './subworkflows/local/extract_feat'

workflow {

    main:

    raw_data_csv = Channel.fromPath(params.input)
    triad_inf_dir = file("$workflow.outputDir/triad/inference").toUriString()
    pmhc_inf_dir = file("$workflow.outputDir/pmhc/inference").toUriString()

    // CLEANING + NEG GENERATION

    CLEAN_CRESTA(raw_data_csv)
    triad_cleaned = CLEAN_CRESTA.out.triad_parquet
    pmhc_cleaned = CLEAN_CRESTA.out.pmhc_parquet

    EXCLUDE_AF3_TRAINING_DATA(triad_cleaned, pmhc_cleaned)
    triad_cleaned = EXCLUDE_AF3_TRAINING_DATA.out.annot_triad_parquet
    triad_excluded = EXCLUDE_AF3_TRAINING_DATA.out.excluded_triad_parquet
    pmhc_cleaned = EXCLUDE_AF3_TRAINING_DATA.out.remaining_pmhc_parquet

    GEN_NEGATIVES(triad_cleaned, triad_cleaned.toList())
    triad_negatives = GEN_NEGATIVES.out.triad_negatives

    // MSA (applies to both unwindowed and windowed sets since peptides aren't aligned)

    MSA_FROM_TRIAD_PARQUET(triad_negatives)
    triad_msa_done_token = MSA_FROM_TRIAD_PARQUET.out.new_msa_list

    // UNWINDOWED (raw peptide) SET

    UNBATCHED_INFERENCE_FROM_UNWIN_TRIAD_PARQUET(triad_negatives, triad_inf_dir, triad_msa_done_token)
    triad_unwin_meta_inf = UNBATCHED_INFERENCE_FROM_UNWIN_TRIAD_PARQUET.out.new_meta_inf
    triad_unwin_neg = DEPEND_TRIAD_ON_UNWIN_INFERENCE(triad_negatives, triad_unwin_meta_inf.toList())

    UNBATCHED_INFERENCE_FROM_UNWIN_PMHC_PARQUET(pmhc_cleaned, pmhc_inf_dir, triad_msa_done_token)
    pmhc_unwin_meta_inf = UNBATCHED_INFERENCE_FROM_UNWIN_PMHC_PARQUET.out.new_meta_inf
    pmhc_unwin = DEPEND_PMHC_ON_UNWIN_INFERENCE(pmhc_cleaned, pmhc_unwin_meta_inf.toList())

    // CONFIDENCE FEATURES
    triad_unwin_conf = EXTRACT_UNWIN_TRIAD_CONF_FEAT(triad_unwin_neg, triad_inf_dir, Channel.value("af3"))
    triad_unwin_tcrdock = TCRDOCK_GEOM_FROM_INFERENCE(triad_unwin_neg, triad_inf_dir, Channel.value("af3"))

    // WINDOWED SET

    // WINDOWING
    PERFORM_WINDOW(triad_negatives)
    triad_neg_w = PERFORM_WINDOW.out.triad_df
    pmhc_neg_w = PERFORM_WINDOW.out.pmhc_df

    // INFERENCE

    UNBATCHED_INFERENCE_FROM_TRIAD_PARQUET(triad_neg_w, triad_inf_dir, triad_msa_done_token)
    triad_meta_inf = UNBATCHED_INFERENCE_FROM_TRIAD_PARQUET.out.new_meta_inf
    triad_win_neg = DEPEND_TRIAD_ON_INFERENCE(triad_neg_w, triad_meta_inf.toList())

    UNBATCHED_INFERENCE_FROM_PMHC_PARQUET(pmhc_neg_w, pmhc_inf_dir, triad_msa_done_token)
    pmhc_meta_inf = UNBATCHED_INFERENCE_FROM_PMHC_PARQUET.out.new_meta_inf
    pmhc_win_neg = DEPEND_PMHC_ON_INFERENCE(pmhc_neg_w, pmhc_meta_inf.toList())

    // CONFIDENCE FEATURES
    triad_conf = EXTRACT_TRIAD_CONF_FEAT(triad_win_neg, triad_inf_dir, Channel.value("af3"))

    pmhc_conf = EXTRACT_PMHC_CONF_FEAT(pmhc_win_neg, pmhc_inf_dir, Channel.value("af3"))


    publish:
    triad_cleaned = triad_cleaned
    pmhc_cleaned = pmhc_cleaned
    triad_excluded = triad_excluded

    // unwindowed
    triad_negatives = triad_negatives
    triad_unwin_meta_inf = triad_unwin_meta_inf
    pmhc_unwin_meta_inf = pmhc_unwin_meta_inf
    triad_unwin_conf = triad_unwin_conf
    triad_unwin_tcrdock = triad_unwin_tcrdock
    

    // windowed
    triad_neg_w = triad_neg_w
    pmhc_neg_w = pmhc_neg_w
    incomplete_negative_log = GEN_NEGATIVES.out.discard_df
    triad_meta_inf = triad_meta_inf
    pmhc_meta_inf = pmhc_meta_inf
    triad_conf = triad_conf
    pmhc_conf = pmhc_conf
    
}

output {
    triad_cleaned {
        path "triad/staged"
    }
    pmhc_cleaned {
        path "pmhc/staged"
    }
    triad_excluded {
        path "triad/staged"
    }
    incomplete_negative_log {
        path "triad/staged"
    }
    triad_negatives {
        path "triad/staged"
    }


    triad_unwin_meta_inf {
        path "triad/inference"
    }
    pmhc_unwin_meta_inf {
        path "pmhc/inference"
    }
    triad_unwin_conf {
        path "triad/staged"
    }
    triad_unwin_tcrdock {
        path "triad/staged"
    }

    triad_neg_w {
        path "triad/staged"
    }
    pmhc_neg_w {
        path "pmhc/staged"
    }

    triad_meta_inf {
        path "triad/inference"
    }
    pmhc_meta_inf {
        path "pmhc/inference"
    }
    triad_conf {
        path "triad/staged"
    }
    pmhc_conf {
        path "pmhc/staged"
    }
}
