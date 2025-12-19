nextflow.preview.output = true
include { BOLTZ_FROM_TRIAD_PARQUET_WITH_TEMPLATES } from '../subworkflows/local/boltz'
include { splitParquet } from 'plugin/nf-parquet'


workflow {

    main:

    input_pq = Channel.fromPath(params.input)

    triad_boltz_dir = file("$workflow.outputDir/triad/predictions_with_templates").toUriString()
    // pmhc_boltz_dir = file("$workflow.outputDir/pmhc/predictions_with_templates").toUriString()

    BOLTZ_FROM_TRIAD_PARQUET_WITH_TEMPLATES(input_pq, triad_boltz_dir)
    triad_all_meta_boltz = BOLTZ_FROM_TRIAD_PARQUET_WITH_TEMPLATES.out.new_meta_inf


    publish:

    triad_all_meta_boltz = triad_all_meta_boltz

}

output {

    triad_all_meta_boltz {
        path "triad/predictions_with_templates"
    }
    

}
