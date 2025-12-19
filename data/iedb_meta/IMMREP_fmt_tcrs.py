import iedb_immrep.utils
import pandas as pd

addtl_xray = pd.read_parquet(
    "addtl_xray.parquet"
)

iedb_immrep.utils.create_leaders_constants("human")

addtl_xray["CDRA3"] = addtl_xray["tcr_1_cdr_3"].str.slice(start=1, stop=-1)
addtl_xray["CDRB3"] = addtl_xray["tcr_2_cdr_3"].str.slice(start=1, stop=-1)
addtl_xray["TRAV"] = addtl_xray["tcr_1_v_gene"]
addtl_xray["TRAJ"] = addtl_xray["tcr_1_j_gene"]
addtl_xray["TRBV"] = addtl_xray["tcr_2_v_gene"]
addtl_xray["TRBJ"] = addtl_xray["tcr_2_j_gene"]
addtl_xray["TCR_name"] = range(addtl_xray.shape[0])

# addtl_xray.rename(
#     columns={
#         "tcr_1_cdr_3": "CDRA3",
#         "tcr_2_cdr_3": "CDRB3",
#         "tcr_1_v_gene": "TRAV",
#         "tcr_1_j_gene": "TRAJ",
#         "tcr_2_v_gene": "TRBV",
#         "tcr_2_j_gene": "TRBJ",
#     },
#     inplace=True,
# )

# addtl_xray["TCR_name"] = addtl_xray["pdb_id"]

thimble_outputs = iedb_immrep.utils.run_thimble(addtl_xray)
addtl_xray["tcr_1_seq"] = addtl_xray["TCR_name"].map(thimble_outputs["A"][0])
addtl_xray["tcr_2_seq"] = addtl_xray["TCR_name"].map(thimble_outputs["B"][0])

addtl_xray = addtl_xray.drop(
    ["CDRA3", "CDRB3", "TRAV", "TRAJ", "TRBV", "TRBJ", "TCR_name", "tcr_1_v_gene", "tcr_2_v_gene", "tcr_1_j_gene", "tcr_2_j_gene"], axis=1
)

addtl_xray.to_parquet("addtl_xray_corrected.parquet", index=False)
