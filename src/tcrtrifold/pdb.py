import requests
from Bio import SeqIO
from io import StringIO
from mdaf3.FeatureExtraction import split_apply_combine
import polars as pl
from tcr_format_parsers.common.MHCCodeConverter import (
    HLASequenceDBConverter,
    H2SequenceDictConverter,
)
import warnings
from pathlib import Path
import MDAnalysis as mda



def get_true_mda_universe(pdb_id, root_path):
    # Favor PDB since it doesn't have multiple residue with same ID issue
    if (root_path / (pdb_id + ".pdb")).exists():
        suffix = ".pdb"

    else:
        suffix = ".cif"

    return mda.Universe((root_path / (pdb_id + suffix)).as_posix())
