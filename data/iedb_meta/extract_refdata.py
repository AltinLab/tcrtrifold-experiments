import polars as pl
import sqlite3

DB_PATH = "/tgen_labs/altin/alphafold3/IEDB/2025-04-15/iedb_public.db"

TRIAD_QUERY_STR = """
SELECT curated_receptor.distinct_receptor_id as receptor_id,
    reference.reference_id,
    assay_type.assay_type
FROM tcell
    JOIN tcell_receptor on tcell.tcell_id = tcell_receptor.tcell_id
    JOIN curated_receptor on tcell_receptor.curated_receptor_id = curated_receptor.curated_receptor_id
    JOIN reference on tcell.reference_id = reference.reference_id
    JOIN assay_type on assay_type.assay_type_id = tcell.as_type_id;
"""

conn = sqlite3.connect(DB_PATH)
iedb_metadat = (
    pl.read_database(
        connection=conn,
        query=TRIAD_QUERY_STR,
        infer_schema_length=1000000,
    )
    .cast({"receptor_id": pl.String, "reference_id": pl.String})
    .select("receptor_id", "reference_id", "assay_type")
    .unique()
)

iedb_metadat.write_parquet("receptor_reference.parquet")
