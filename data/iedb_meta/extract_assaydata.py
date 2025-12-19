import polars as pl
import sqlite3
from tcrtrifold.utils import generate_job_name
from tcrtrifold.cleaning import infer_hla_chain

DB_PATH = "/tgen_labs/altin/alphafold3/IEDB/2025-04-15/iedb_public.db"

TRIAD_QUERY_STR = """
SELECT distinct_receptor.distinct_receptor_id as receptor_id, 
    tcell.tcell_id, mhc_allele_name, mhc_allele_restriction.class as mhc_class,
    assay_type.assay_type, reference.reference_id, epitope.linear_peptide_seq
FROM tcell
    JOIN tcell_receptor on tcell.tcell_id = tcell_receptor.tcell_id
    JOIN assay_type on assay_type.assay_type_id = tcell.as_type_id
    JOIN curated_receptor on tcell_receptor.curated_receptor_id = curated_receptor.curated_receptor_id
    JOIN distinct_receptor on curated_receptor.distinct_receptor_id = distinct_receptor.distinct_receptor_id
    JOIN mhc_allele_restriction on mhc_allele_restriction.mhc_allele_restriction_id = tcell.mhc_allele_restriction_id
    JOIN reference on tcell.reference_id = reference.reference_id
    JOIN curated_epitope on tcell.curated_epitope_id = curated_epitope.curated_epitope_id
    JOIN epitope_object on curated_epitope.e_object_id = epitope_object.object_id
    JOIN epitope on epitope.epitope_id = epitope_object.epitope_id;
"""

conn = sqlite3.connect(DB_PATH)
iedb_metadat = (
    pl.read_database(
        connection=conn,
        query=TRIAD_QUERY_STR,
        infer_schema_length=1000000,
    )
    .cast({"receptor_id": pl.String, "reference_id": pl.String})
    .rename({"reference_id": "references"})
    # .select("receptor_id", "assay_type")
    .unique()
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
)

iedb_II = (
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
                pl.col("split_parts").list.get(1, null_on_oob=True).alias("mhc_2_name"),
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
).select(pl.exclude("mhc_struct", "split_parts", "mhc_name_inferred"))

print(iedb_I.columns, iedb_II.columns)

iedb_metadat = pl.concat([iedb_I, iedb_II])

iedb_metadat.write_parquet("assay_type.parquet")
