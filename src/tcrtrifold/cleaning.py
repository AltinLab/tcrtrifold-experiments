from .mhc import (
    shorten_to_fullname,
    is_fullname,
    DQA_FOR,
    DPA_FOR,
)
import warnings
from .utils import (
    SOURCE_RENAME_DICT,
    SOURCE_ANTIGEN_COLS,
    TCRDIST_COLS,
)
import random
from tqdm import tqdm
import polars as pl
import MDAnalysis as mda
from mdaf3.FeatureExtraction import *
from mdaf3.AF3OutputParser import *
import Bio.Align


def infer_hla_chain(mhc_1_name, mhc_2_name):

    nullchains = {
        "mhc_1_name": None,
        "mhc_2_name": None,
        "mhc_name_inferred": None,
    }

    if mhc_2_name.startswith("DRB"):
        fullname = shorten_to_fullname(mhc_2_name)
        if is_fullname(fullname):
            return {
                # use 0102 sicne sequence in uniprot
                # mutation outside top region
                "mhc_1_name": "DRA*01:02",
                "mhc_2_name": fullname,
                # only one option- we don't count this as inferred
                "mhc_name_inferred": "neither",
            }

        else:
            warnings.warn(f"Could not find fullname for {mhc_2_name}")
            return nullchains
    elif mhc_2_name.startswith("DQB"):
        fullname = shorten_to_fullname(mhc_2_name)
        if is_fullname(fullname):

            b_chain = fullname
            a_chain = (
                shorten_to_fullname(mhc_1_name) if mhc_1_name is not None else None
            )
            if mhc_1_name is None or not is_fullname(a_chain):
                if b_chain in DQA_FOR:
                    a_chain = DQA_FOR[b_chain]
                    return {
                        "mhc_1_name": a_chain,
                        "mhc_2_name": b_chain,
                        "mhc_name_inferred": "chain_1",
                    }
                else:
                    warnings.warn(f"Could not find DQA chain for {b_chain}")
                    return nullchains
            else:
                return {
                    "mhc_1_name": a_chain,
                    "mhc_2_name": b_chain,
                    "mhc_name_inferred": "neither",
                }
        else:
            return nullchains

    elif mhc_2_name.startswith("DPB"):
        fullname = shorten_to_fullname(mhc_2_name)
        if is_fullname(fullname):
            b_chain = fullname
            a_chain = (
                shorten_to_fullname(mhc_1_name) if mhc_1_name is not None else None
            )
            if mhc_1_name is None or not is_fullname(a_chain):

                if b_chain in DPA_FOR:
                    a_chain = DPA_FOR[b_chain]
                    return {
                        "mhc_1_name": a_chain,
                        "mhc_2_name": b_chain,
                        "mhc_name_inferred": "chain_1",
                    }
                else:
                    warnings.warn(f"Could not find DPA chain for {b_chain}")
                    return nullchains
            else:
                return {
                    "mhc_1_name": a_chain,
                    "mhc_2_name": b_chain,
                    "mhc_name_inferred": "neither",
                }
        else:
            return nullchains
    else:
        return nullchains
