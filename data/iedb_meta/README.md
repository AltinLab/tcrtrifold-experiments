
Assay type data was extracted from IEDB using a local clone of IEDB (since the HTTP API doesn't expose the underlying schema) converted to SQLite (since persistent DB servers are non-trivial on TGen HPC).

To replicate this process, use this page to export IEDB to MySQL: https://www.iedb.org/database_export_v3.php

And then use this tool to convert the MySQL dump to SQLite: https://github.com/mysql2sqlite/mysql2sqlite

See `extract_assaydata.py` for the query we used.