def tcrdock_format_cif(row, inference_path, output_path):
    if Path(output_path / (row["job_name"] + ".pdb")).exists():
        return pl.DataFrame([row])

    af3_output = AF3Output(inference_path / row["job_name"])

    pred_u = af3_output.get_mda_universe()

    # maybe mhc 2 seq is not included since it's implied B2m, so check
    if row["mhc_class"] == "II":
        # mhc1, mhc2, pep, tcr1, tcr2
        pred_segids = ["B", "C", "A", "D", "E"]
        rename_segids = ["A", "B", "C", "D", "E"]
    # TCRdock will remove B2M anyways
    else:
        pred_segids = ["B", "A", "D", "E"]
        # for consistency with format already in tcrdock repo
        rename_segids = ["A", "B", "C", "D"]

    chain_us = []

    for pred_segsel, rename_segid in zip(pred_segids, rename_segids):

        pred_sel = pred_u.select_atoms(f"segid {pred_segsel}").atoms

        chain_u = mda.Merge(pred_sel)
        chain_u.segments.segids = rename_segid
        chain_u.atoms.chainIDs = [rename_segid] * len(chain_u.atoms)
        chain_us.append(chain_u.atoms)

    new_u = mda.Merge(*chain_us)

    with mda.Writer(output_path / (row["job_name"] + ".pdb")) as W:

        # u_new = mda.Universe.empty(
        #     n_atoms, n_segments=n_segments, n_residues=n_residues
        # )

        # ordered_chains = sum(chain_sels)

        # # for attr in ["name", "type", "resname"]:
        # #     u_new.add_TopologyAttr("name", ordered_chains.residues.names)

        # # choose first altloc if mutliple present
        # W.write(ordered_chains)

        W.write(new_u.atoms)

    # noop
    return row
