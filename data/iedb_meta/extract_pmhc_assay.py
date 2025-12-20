import polars as pl
import sqlite3
from tcrtrifold.utils import generate_job_name
from tcrtrifold.cleaning import infer_hla_chain

DB_PATH = "/tgen_labs/altin/alphafold3/IEDB/2025-04-15/iedb_public.db"

TRIAD_QUERY_STR = """
SELECT mhc_allele_name,
    linear_peptide_seq,
    assay_type.*,
    mhc_bind.mhc_bind_id,
    mhc_bind.reference_id,
    mhc_bind.as_char_value,
    mhc_bind.as_num_value,
    mhc_bind.as_location,
    mhc_bind.as_type_id,
    mhc_bind.as_comments,
    mhc_allele_restriction.class as mhc_class
FROM mhc_bind
    JOIN assay_type on assay_type.assay_type_id = mhc_bind.as_type_id
    JOIN mhc_allele_restriction on mhc_allele_restriction.mhc_allele_restriction_id = mhc_bind.mhc_allele_restriction_id
    JOIN curated_epitope on curated_epitope.curated_epitope_id = mhc_bind.curated_epitope_id
    JOIN epitope_object on curated_epitope.e_object_id = epitope_object.object_id
    JOIN epitope on epitope.epitope_id = epitope_object.epitope_id;
"""

conn = sqlite3.connect(DB_PATH)
iedb_metadat = (
    (
        pl.read_database(
            connection=conn,
            query=TRIAD_QUERY_STR,
            infer_schema_length=1000000,
        )
        # .select("receptor_id", "assay_type")
        .unique()
    )
    .filter(pl.col("response").str.contains("KD"))
    .filter(~pl.col("mhc_allele_name").str.contains("mutant"))
    .filter(~pl.col("mhc_allele_name").is_null())
)
iedb_I = (
    iedb_metadat.filter(pl.col("mhc_class") == "I")
    .rename({"mhc_allele_name": "mhc_1_name", "linear_peptide_seq": "peptide"})
    .with_columns(
        pl.col("mhc_1_name").str.split_exact("HLA-", 1).alias("split_parts"),
    )
    .select(pl.exclude("mhc_1_name"))
    .unnest("split_parts")
    .rename(
        {
            "field_0": "tmp",
            "field_1": "mhc_1_name",
        }
    )
    .select(pl.exclude("tmp", "mhc_struct", "split_parts", "mhc_name_inferred"))
    .with_columns(pl.lit("B2M").alias("mhc_2_name"))
).filter(~pl.col("mhc_1_name").is_null())

iedb_II = (
    (
        iedb_metadat.filter(pl.col("mhc_class") == "II")
        .rename({"mhc_allele_name": "mhc_1_name", "linear_peptide_seq": "peptide"})
        .with_columns(pl.col("mhc_1_name").str.split("/").alias("split_parts"))
        .with_columns(
            pl.when(pl.col("split_parts").list.len() == 2)
            .then(
                pl.struct(
                    pl.col("split_parts")
                    .list.get(0, null_on_oob=True)
                    .str.slice(4)
                    .alias("mhc_1_name"),
                    pl.col("split_parts")
                    .list.get(1, null_on_oob=True)
                    .alias("mhc_2_name"),
                )
            )
            .otherwise(
                pl.struct(
                    pl.lit(None).alias("mhc_1_name"),
                    pl.col("split_parts").list.get(0).str.slice(4).alias("mhc_2_name"),
                )
            )
            .alias("mhc_struct")
        )
        .select(pl.exclude("mhc_1_name"))
        .with_columns(
            pl.col("mhc_struct")
            .map_elements(
                lambda x: infer_hla_chain(x["mhc_1_name"], x["mhc_2_name"]),
                return_dtype=pl.Struct(
                    {
                        "mhc_1_name": pl.String,
                        "mhc_2_name": pl.String,
                        "mhc_name_inferred": pl.String,
                    }
                ),
            )
            .alias("chains")
        )
        .unnest("chains")
    )
    .select(pl.exclude("mhc_struct", "split_parts", "mhc_name_inferred"))
    .filter(~pl.col("mhc_2_name").is_null())
)

print(iedb_I.columns, iedb_II.columns)

iedb_metadat = pl.concat([iedb_I, iedb_II])

iedb_metadat.write_csv("pmhc_assays.tsv", separator="\t")
