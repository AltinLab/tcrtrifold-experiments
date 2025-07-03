def find_matching(row, blast_result, thresh=95):
    """
    A function to find matching PDBs for a given peptide and MHC
    sequence. If the peptide, MHC heavy (I) or MHC alpha and beta (II)
    all match a PDB entry with a per_identity greater than the
    threshold, then the PDB entry is considered a match.
    """
    pep_match = (
        blast_result.filter(
            pl.col("query_id") == row["peptide"],
            pl.col("per_identity") >= thresh,
        )
        .select("pdb", "per_identity")
        .rename({"per_identity": "peptide_per_identity"})
        .unique()
    )
    mhc_1_match = (
        blast_result.filter(
            pl.col("query_id") == row["mhc_1_seq"],
            pl.col("per_identity") >= thresh,
        )
        .select("pdb", "per_identity")
        .rename({"per_identity": "mhc_1_per_identity"})
        .unique()
    )

    if row["mhc_class"] == "II":
        mhc_2_match = (
            blast_result.filter(
                pl.col("query_id") == row["mhc_2_seq"],
                pl.col("per_identity") >= thresh,
            )
            .select("pdb", "per_identity")
            .rename({"per_identity": "mhc_2_per_identity"})
            .unique()
        )

        pmhc_match_ids = (
            pep_match.join(mhc_1_match, on="pdb")
            .join(mhc_2_match, on="pdb")
            .sort(
                by=[
                    "peptide_per_identity",
                    "mhc_1_per_identity",
                    "mhc_2_per_identity",
                ],
                descending=True,
            )
        )

    else:
        pmhc_match_ids = pep_match.join(mhc_1_match, on="pdb").sort(
            by=[
                "peptide_per_identity",
                "mhc_1_per_identity",
            ],
            descending=True,
        )

    pmhc_match = pmhc_match_ids.height > 0

    new_row = row.copy()

    new_row["pmhc_match"] = pmhc_match
    new_row["pmhc_match_id"] = (
        pmhc_match_ids.select("pdb").to_series().to_list()
        if pmhc_match
        else None
    )

    if len(pmhc_match_ids.select("pdb").to_series().to_list()) > 1:
        print(row["mhc_1_name"])

    new_row["pmhc_match_mhc_1_per_identity"] = (
        pmhc_match_ids[0].select("mhc_1_per_identity").item()
        if pmhc_match
        else None
    )
    new_row["pmhc_match_mhc_2_per_identity"] = (
        pmhc_match_ids[0].select("mhc_2_per_identity").item()
        if pmhc_match and row["mhc_class"] == "II"
        else None
    )
    new_row["pmhc_match_peptide_per_identity"] = (
        pmhc_match_ids[0].select("peptide_per_identity").item()
        if pmhc_match
        else None
    )
    return pl.DataFrame([new_row])
