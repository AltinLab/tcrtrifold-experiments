import polars as pl
import sqlite3

DB_PATH = "/tgen_labs/altin/alphafold3/IEDB/2025-04-15/iedb_public.db"

TRIAD_QUERY_STR = """
SELECT article.reference_id,
    article.article_date,
    article.medline_date
FROM article
"""

conn = sqlite3.connect(DB_PATH)
iedb_metadat = pl.read_database(
    connection=conn,
    query=TRIAD_QUERY_STR,
    infer_schema_length=1000000,
).cast({"reference_id": pl.String})

iedb_metadat.write_parquet("reference_date.parquet")
