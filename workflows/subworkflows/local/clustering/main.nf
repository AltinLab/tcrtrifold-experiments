include { SPLIT_TRIAD_INTO_CHAINS as SPLIT_QUERY; 
    SPLIT_TRIAD_INTO_CHAINS as SPLIT_VALIDATION;
    PARQUET_TO_FASTA as PEPTIDE_TO_FASTA_QUERY;
    PARQUET_TO_FASTA as MHC_1_TO_FASTA_QUERY;
    PARQUET_TO_FASTA as MHC_2_TO_FASTA_QUERY;
    PARQUET_TO_FASTA as TCR_1_TO_FASTA_QUERY;
    PARQUET_TO_FASTA as TCR_2_TO_FASTA_QUERY;
    PARQUET_TO_FASTA as PEPTIDE_TO_FASTA_VALIDATION;
    PARQUET_TO_FASTA as MHC_1_TO_FASTA_VALIDATION;
    PARQUET_TO_FASTA as MHC_2_TO_FASTA_VALIDATION;
    PARQUET_TO_FASTA as TCR_1_TO_FASTA_VALIDATION;
    PARQUET_TO_FASTA as TCR_2_TO_FASTA_VALIDATION;
    MMSEQS_QUERY_TARGET_IDENT as MMSEQS_PEPTIDE;
    MMSEQS_QUERY_TARGET_IDENT as MMSEQS_MHC_1;
    MMSEQS_QUERY_TARGET_IDENT as MMSEQS_MHC_2;
    MMSEQS_QUERY_TARGET_IDENT as MMSEQS_TCR_1;
    MMSEQS_QUERY_TARGET_IDENT as MMSEQS_TCR_2;
    FILTER_ANNOTATE_TRIAD_FROM_IDENT;
    FILTER_ANNOTATE_TRIAD_FROM_PDB; 
    MMSEQS_QUERY_TARGET_IDENT_SHORT;
    BLAST_PARQUET_WEBSERVER as BLAST_MHC_1;
    BLAST_PARQUET_WEBSERVER as BLAST_MHC_2;
    BLAST_PARQUET_WEBSERVER as BLAST_TCR_1;
    BLAST_PARQUET_WEBSERVER as BLAST_TCR_2;
} from '../../../modules/local/clustering'

workflow EXCLUDE_VALIDATION_TRIADS {
    /*
    Perform two tasks:
    - annotate the antigens that are present in the validation set in the query set
    - REMOVE the triads from the query set that are presesnt in the validation set
    */

    take:
    query_triad_parquet
    curr_pmhc_parquet
    validation_triad_parquet

    main:

    SPLIT_QUERY(query_triad_parquet)
    SPLIT_VALIDATION(validation_triad_parquet)

    PEPTIDE_TO_FASTA_QUERY(SPLIT_QUERY.out.peptide_pq)
    MHC_1_TO_FASTA_QUERY(SPLIT_QUERY.out.mhc_1_pq)
    MHC_2_TO_FASTA_QUERY(SPLIT_QUERY.out.mhc_2_pq)
    TCR_1_TO_FASTA_QUERY(SPLIT_QUERY.out.tcr_1_pq)
    TCR_2_TO_FASTA_QUERY(SPLIT_QUERY.out.tcr_2_pq)

    PEPTIDE_TO_FASTA_VALIDATION(SPLIT_VALIDATION.out.peptide_pq)
    MHC_1_TO_FASTA_VALIDATION(SPLIT_VALIDATION.out.mhc_1_pq)
    MHC_2_TO_FASTA_VALIDATION(SPLIT_VALIDATION.out.mhc_2_pq)
    TCR_1_TO_FASTA_VALIDATION(SPLIT_VALIDATION.out.tcr_1_pq)
    TCR_2_TO_FASTA_VALIDATION(SPLIT_VALIDATION.out.tcr_2_pq)

    MMSEQS_QUERY_TARGET_IDENT_SHORT(PEPTIDE_TO_FASTA_QUERY.out.fasta, 
        PEPTIDE_TO_FASTA_VALIDATION.out.fasta)

    MMSEQS_MHC_1(MHC_1_TO_FASTA_QUERY.out.fasta,
        MHC_1_TO_FASTA_VALIDATION.out.fasta)

    MMSEQS_MHC_2(MHC_2_TO_FASTA_QUERY.out.fasta,
        MHC_2_TO_FASTA_VALIDATION.out.fasta)

    MMSEQS_TCR_1(TCR_1_TO_FASTA_QUERY.out.fasta,
        TCR_1_TO_FASTA_VALIDATION.out.fasta)

    MMSEQS_TCR_2(TCR_2_TO_FASTA_QUERY.out.fasta,
        TCR_2_TO_FASTA_VALIDATION.out.fasta)

    FILTER_ANNOTATE_TRIAD_FROM_IDENT(
        SPLIT_QUERY.out.triad_pq,
        curr_pmhc_parquet,
        SPLIT_VALIDATION.out.triad_pq,
        MMSEQS_QUERY_TARGET_IDENT_SHORT.out.ident_tsv,
        MMSEQS_MHC_1.out.ident_tsv,
        MMSEQS_MHC_2.out.ident_tsv,
        MMSEQS_TCR_1.out.ident_tsv,
        MMSEQS_TCR_2.out.ident_tsv
    )


    emit:
    annot_triad_parquet = FILTER_ANNOTATE_TRIAD_FROM_IDENT.out.annot_triad_parquet
    excluded_triad_parquet = FILTER_ANNOTATE_TRIAD_FROM_IDENT.out.excluded_triad_parquet
    remaining_pmhc_parquet = FILTER_ANNOTATE_TRIAD_FROM_IDENT.out.remaining_pmhc_parquet

}

workflow EXCLUDE_AF3_TRAINING_DATA {
    take:
    query_triad_parquet
    curr_pmhc_parquet

    main:

    SPLIT_QUERY(query_triad_parquet)

    BLAST_MHC_1(SPLIT_QUERY.out.mhc_1_pq)
    BLAST_MHC_2(SPLIT_QUERY.out.mhc_2_pq)
    BLAST_TCR_1(SPLIT_QUERY.out.tcr_1_pq)
    BLAST_TCR_2(SPLIT_QUERY.out.tcr_2_pq)

    FILTER_ANNOTATE_TRIAD_FROM_PDB(
        SPLIT_QUERY.out.triad_pq,
        curr_pmhc_parquet,
        BLAST_MHC_1.out.blast_tsv,
        BLAST_MHC_2.out.blast_tsv,
        BLAST_TCR_1.out.blast_tsv,
        BLAST_TCR_2.out.blast_tsv
    )

    emit:
    annot_triad_parquet = FILTER_ANNOTATE_TRIAD_FROM_PDB.out.annot_triad_parquet
    excluded_triad_parquet = FILTER_ANNOTATE_TRIAD_FROM_PDB.out.excluded_triad_parquet
    remaining_pmhc_parquet = FILTER_ANNOTATE_TRIAD_FROM_PDB.out.remaining_pmhc_parquet
}