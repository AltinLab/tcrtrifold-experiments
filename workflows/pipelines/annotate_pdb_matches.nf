nextflow.preview.output = true

include { SPLIT_TRIAD_INTO_CHAINS as SPLIT_QUERY; 
    ANNOTATE_TRIAD_FROM_PDB; 
    MMSEQS_QUERY_TARGET_IDENT_SHORT;
    BLAST_PARQUET_WEBSERVER as BLAST_MHC_1;
    BLAST_PARQUET_WEBSERVER as BLAST_MHC_2;
    BLAST_PARQUET_WEBSERVER as BLAST_TCR_1;
    BLAST_PARQUET_WEBSERVER as BLAST_TCR_2;
} from '../modules/local/clustering'


workflow {

    main:

    SPLIT_QUERY(Channel.fromPath(params.input))

    BLAST_MHC_1(SPLIT_QUERY.out.mhc_1_pq)
    BLAST_MHC_2(SPLIT_QUERY.out.mhc_2_pq)
    BLAST_TCR_1(SPLIT_QUERY.out.tcr_1_pq)
    BLAST_TCR_2(SPLIT_QUERY.out.tcr_2_pq)

    ANNOTATE_TRIAD_FROM_PDB(
        SPLIT_QUERY.out.triad_pq,
        BLAST_MHC_1.out.blast_tsv,
        BLAST_MHC_2.out.blast_tsv,
        BLAST_TCR_1.out.blast_tsv,
        BLAST_TCR_2.out.blast_tsv
    )

    publish:
    annot_triad_parquet = ANNOTATE_TRIAD_FROM_PDB.out.annot_triad_parquet
}

output {

    annot_triad_parquet { path "triad" }
}