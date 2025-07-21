from torch_geometric.data import Data
from torch_geometric.transforms import ToUndirected
import torch_geometric.transforms as T
from torch_geometric.utils import remove_isolated_nodes
import h5py
from MDAnalysis.analysis.distances import distance_array, self_distance_array
from scipy.spatial.distance import squareform
import torch
import numpy as np


def mda_triad_to_graph(u, contact_prob, mhc_class, cognate, id):

    focal_atm = u.select_atoms("segid A or around 5 segid A")
    focal_res = focal_atm.residues
    focal_calpha = focal_res.atoms.select_atoms("name CA")

    node_positions = torch.from_numpy(focal_calpha.positions)
    node_confidences = torch.tensor([res.atoms.tempfactors.mean() for res in focal_res])

    edge_index_list = []

    for i, res_1 in enumerate(focal_res[:-1]):
        for j, res_2 in zip(range(i + 1, len(focal_res)), focal_res[i + 1 :]):
            min_d = distance_array(res_1.atoms.positions, res_2.atoms.positions).min()
            if min_d < 5:
                edge_index_list.append([i, j])

    edge_index = np.array(edge_index_list)
    # use calpha to calculate edge dist
    # can result in edge >5 angstrom
    edge_distances = squareform(self_distance_array(focal_calpha.positions))[
        edge_index[:, 0], edge_index[:, 1]
    ]

    row_residx = focal_res[edge_index[:, 0]].resindices
    col_residx = focal_res[edge_index[:, 1]].resindices
    edge_contact_prob = contact_prob[row_residx, col_residx]

    x = torch.from_numpy(np.stack([node_confidences]).T)
    edge_attr = torch.from_numpy(np.stack([edge_distances, edge_contact_prob]).T)

    edge_index = torch.from_numpy(edge_index.T)

    # edge_index, edge_attr, node_boolmask = remove_isolated_nodes(
    #     edge_index, edge_attr, num_nodes=x.shape[0]
    # )

    # x = x[node_boolmask]
    # node_positions = node_positions[node_boolmask]

    segids = list(str(ch) for ch in focal_res.segids.astype(np.str_))

    graph = Data(
        x=x,
        edge_index=edge_index,
        y=torch.tensor([cognate]),
        edge_attr=edge_attr,
        pos=node_positions,
        id=id,
        segids=segids,
    )

    return ToUndirected()(graph)
