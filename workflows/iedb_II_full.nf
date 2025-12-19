/*
Parameters:
- input: Path to the IEDB II raw data CSV file
*/

// in current version, new output syntax is in preview
nextflow.preview.output = true


include { MSA_WORKFLOW } from './subworkflows/tgen/af3'
include { SEQ_LIST_TO_FASTA } from './modules/tgen/af3'
include { CLEAN_IEDB_II } from './subworkflows/local/cleaning'
include { EXCLUDE_VALIDATION_TRIADS } from './subworkflows/local/clustering'
include { GEN_NEGATIVES; PERFORM_WINDOW } from './subworkflows/local/neg'
include { MSA_FROM_TRIAD_PARQUET; 
        UNBATCHED_INFERENCE_FROM_TRIAD_PARQUET; 
        UNBATCHED_INFERENCE_FROM_PMHC_PARQUET;
        NOOP_DEP as DEPEND_TRIAD_ON_INFERENCE; 
        NOOP_DEP as DEPEND_PMHC_ON_INFERENCE  } from './subworkflows/local/af3_adapter'
include { EXTRACT_TRIAD_CONF_FEAT; TCRDOCK_GEOM_FROM_INFERENCE; EXTRACT_PMHC_CONF_FEAT } from './subworkflows/local/extract_feat'

workflow {

    main:

    raw_data_csv = Channel.fromPath(params.input)
    triad_inf_dir = file("$workflow.outputDir/triad/inference").toUriString()
    pmhc_inf_dir = file("$workflow.outputDir/pmhc/inference").toUriString()


    CLEAN_IEDB_II(raw_data_csv)
    triad_cleaned = CLEAN_IEDB_II.out.triad_parquet
    pmhc_cleaned = CLEAN_IEDB_II.out.pmhc_parquet

    GEN_NEGATIVES(triad_cleaned, triad_cleaned.toList())
    triad_negatives = GEN_NEGATIVES.out.triad_negatives


    MSA_FROM_TRIAD_PARQUET(triad_negatives)
    triad_msa_done_token = MSA_FROM_TRIAD_PARQUET.out.new_msa_list

    UNBATCHED_INFERENCE_FROM_TRIAD_PARQUET(triad_negatives, triad_inf_dir, triad_msa_done_token)
    triad_meta_inf = UNBATCHED_INFERENCE_FROM_TRIAD_PARQUET.out.new_meta_inf
    triad_negatives = DEPEND_TRIAD_ON_INFERENCE(triad_negatives, triad_meta_inf.toList())

    UNBATCHED_INFERENCE_FROM_PMHC_PARQUET(pmhc_cleaned, pmhc_inf_dir, triad_msa_done_token)
    pmhc_meta_inf = UNBATCHED_INFERENCE_FROM_PMHC_PARQUET.out.new_meta_inf
    pmhc_cleaned = DEPEND_PMHC_ON_INFERENCE(pmhc_cleaned, pmhc_meta_inf.toList())

    triad_tcrdock = TCRDOCK_GEOM_FROM_INFERENCE(triad_negatives, triad_inf_dir, Channel.value("af3"))

    publish:
    triad_cleaned = triad_cleaned
    pmhc_cleaned = pmhc_cleaned

    incomplete_negative_log = GEN_NEGATIVES.out.discard_df

    triad_negatives = triad_negatives
    triad_meta_inf = triad_meta_inf
    pmhc_meta_inf = pmhc_meta_inf
    triad_tcrdock = triad_tcrdock
    
}

output {
    triad_cleaned {
        path "triad/staged"
    }
    pmhc_cleaned {
        path "pmhc/staged"
    }
    triad_negatives {
        path "triad/staged"
    }
    incomplete_negative_log {
        path "triad/staged"
    }
    triad_meta_inf {
        path "triad/inference"
    }
    pmhc_meta_inf {
        path "pmhc/inference"
    }
    triad_tcrdock {
        path "triad/staged"
    }
}
