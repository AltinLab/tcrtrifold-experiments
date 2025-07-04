import anarci
import numpy as np


def extract_tcrdist_cols(row):

    for tcr in [1, 2]:
        regions = tcr_by_imgt_region(
            row[f"tcr_{tcr}_seq"],
            np.arange(len(row[f"tcr_{tcr}_seq"])),
            row[f"tcr_{tcr}_chain"],
            row[f"tcr_{tcr}_species"],
        )
        for region in ["cdr_1", "cdr_2", "cdr_2_5", "cdr_3"]:
            row[f"tcr_{tcr}_{region}"] = "".join(
                np.array(list(row[f"tcr_{tcr}_seq"]))[regions[region]]
            )

        v_gene, j_gene = tcr_v_j_genes(
            row[f"tcr_{tcr}_seq"],
            row[f"tcr_{tcr}_chain"],
            row[f"tcr_{tcr}_species"],
        )
        row[f"tcr_{tcr}_v_gene"] = v_gene
        row[f"tcr_{tcr}_j_gene"] = j_gene

    return row


def annotate_tcr(tcr_seq, resindices, tcr_chain, species, strict=False):
    """Return the resindices of a contiguous substring of resindices
    which aligned to a TCR seq using ANARCI"""
    if len(tcr_seq) != len(resindices):
        raise ValueError(
            "TCR sequence length must match length of " "residue index array"
        )

    if tcr_chain == "alpha":
        # quirk of ANARCI
        # https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabpred/anarci/
        allow = set(["A", "D"])
    elif tcr_chain == "beta":
        allow = set(["B"])
    else:
        raise ValueError(f"Invalid chain '{tcr_chain}'")

    results = anarci.run_anarci(
        [("AAA", tcr_seq)],
        scheme="imgt",
        allowed_species=[species],
        allow=allow,
    )
    # check for None
    if results[1][0] is None:
        raise ValueError(f"No domain found for sequence '{tcr_seq}'")
    elif len(results[1][0]) > 1:
        raise ValueError("Multiple domains found for sequence")

    # validate that the species provided is what anarci found
    if strict:
        if results[2][0][0]["species"] != species:
            raise ValueError(f"Species mismatch: {species} vs {results[1][0][0][1]}")

    numbering = results[1][0][0]

    sub_start = numbering[1]
    sub_stop = numbering[2] + 1

    # this entire substring will be numbered
    # however, there may be gaps in the sequence
    # which are given an IMGT number
    numbered_substring = tcr_seq[sub_start:sub_stop]
    resindices_slice = resindices[sub_start:sub_stop]
    imgt_num = np.zeros((len(range(sub_start, sub_stop)),), dtype=np.int32)
    imgt_tuples = numbering[0]
    j = 0
    for i in range(len(numbered_substring)):
        aa = numbered_substring[i]
        while j < len(imgt_tuples) and imgt_tuples[j][1] != aa:
            j += 1
        if j < len(imgt_tuples):
            imgt_num[i] = imgt_tuples[j][0][0]
            j += 1

    # zero is not an IMGT number, so we use this as a quick error cehck
    n_zeroes = np.count_nonzero(imgt_num == 0)
    if n_zeroes != 0:
        raise ValueError(
            f"0 is not an IMGT number. Numbering failed for TCR '{tcr_seq}'"
        )

    return resindices_slice, imgt_num, (sub_start, sub_stop)


def tcr_by_imgt_region(tcr_seq, tcr_resindices, tcr_chain, tcr_species):
    tcr_indices, imgt_num, _ = annotate_tcr(
        tcr_seq,
        tcr_resindices,
        tcr_chain,
        tcr_species,
    )

    tcr_residx_dict = {}

    tcr_residx_dict["fwr_1"] = tcr_indices[((imgt_num >= 1) & (imgt_num <= 26))]

    tcr_residx_dict["cdr_1"] = tcr_indices[((imgt_num >= 27) & (imgt_num <= 38))]

    tcr_residx_dict["fwr_2"] = tcr_indices[((imgt_num >= 39) & (imgt_num <= 55))]

    tcr_residx_dict["cdr_2"] = tcr_indices[((imgt_num >= 56) & (imgt_num <= 65))]

    tcr_residx_dict["fwr_3"] = tcr_indices[((imgt_num >= 66) & (imgt_num <= 103))]

    tcr_residx_dict["cdr_2_5"] = tcr_indices[((imgt_num >= 81) & (imgt_num <= 86))]

    tcr_residx_dict["cdr_3"] = tcr_indices[((imgt_num >= 104) & (imgt_num <= 118))]

    tcr_residx_dict["fwr_4"] = tcr_indices[((imgt_num >= 119) & (imgt_num <= 129))]

    return tcr_residx_dict


def tcr_v_j_genes(tcr_seq, tcr_chain, species):

    if tcr_chain == "alpha":
        # quirk of ANARCI
        # https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabpred/anarci/
        allow = set(["A", "D"])
    elif tcr_chain == "beta":
        allow = set(["B"])
    else:
        raise ValueError(f"Invalid chain '{tcr_chain}'")

    genes = anarci.run_anarci(
        [("AAA", tcr_seq)],
        scheme="imgt",
        allowed_species=[species],
        allow=allow,
        assign_germline=True,
    )[2][0][0]["germlines"]

    return genes["v_gene"][0][1], genes["j_gene"][0][1]
