from mdaf3.AF3OutputParser import AF3Output
import gzip
import MDAnalysis as mda
from MDAnalysis.lib.util import convert_aa_code
import polars as pl
import numpy as np
from pathlib import Path
import argparse


def res_signed_vol(res):
    resname = convert_aa_code(res.resname)

    if resname != "G":
        N = res.atoms.select_atoms("name N").positions[0]
        CA = res.atoms.select_atoms("name CA").positions[0]
        C = res.atoms.select_atoms("name C").positions[0]
        CB = res.atoms.select_atoms("name CB").positions[0]

        V_alpha = float(np.dot(np.cross(N - CA, C - CA), CB - CA))
    else:
        V_alpha = None

    if resname in "IT":
        sel1 = "OG1" if resname == "T" else "CG1"
        CA = res.atoms.select_atoms("name CA").positions[0]
        CB = res.atoms.select_atoms("name CB").positions[0]
        CG1 = res.atoms.select_atoms(f"name {sel1}").positions[0]
        CG2 = res.atoms.select_atoms("name CG2").positions[0]

        V_beta = float(np.dot(np.cross(CA - CB, CG1 - CB), CG2 - CB))
    else:
        V_beta = None

    def sgn(v):
        if v is None:
            return 0
        if v > 1e-4:
            return 1
        if v < -1e-4:
            return -1
        return 0

    return sgn(V_alpha), sgn(V_beta)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_parquet",
        type=str,
    )
    parser.add_argument(
        "--inf_dir",
        type=str,
    )
    parser.add_argument(
        "--output_dir",
        type=str,
    )
    args = parser.parse_args()

    inf_dir = Path(args.inf_dir)

    pq = pl.read_parquet(args.input_parquet)

    chiral_df = []
    for row in pq.iter_rows(named=True):

        for seed in [1, 2, 3, 4, 5]:

            af3_out = AF3Output(inf_dir / row["job_name"])

            u = af3_out.get_mda_universe(seed=seed)

            df = []
            for res in u.residues:

                sgn_a, sgn_b = res_signed_vol(res)

                df.append(
                    {
                        "res": convert_aa_code(res.resname),
                        "resindex": res.resindex,
                        "segid": res.segid,
                        "V_alpha": sgn_a,
                        "V_beta": sgn_b,
                        "seed": seed,
                    }
                )

            df = pl.DataFrame(df)

            # exp_alpha_pos = df.select(pl.col("V_alpha").sum()).item() >0
            # exp_beta_pos = df.select(pl.col("V_beta").sum()).item() >0

            chiral_D = (
                df.filter(pl.col("V_alpha") < 0).select("resindex", "segid").to_numpy()
            )

            row["chiral_D_segid"] = chiral_D[:, 1].tolist()
            row["chiral_D_resindex"] = chiral_D[:, 0].tolist()

            chiral_IT_D = df.filter(pl.col("V_beta") < 0).select("resindex", "segid")

            row["chiral_IT_D_segid"] = chiral_D[:, 1].tolist()
            row["chiral_IT_D_resindex"] = chiral_D[:, 0].tolist()

            chiral_df.append(row)

    chiral_df = pl.DataFrame(chiral_df)

    # if chiral_atoms:
    #     print("hello world")
