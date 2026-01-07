from nvalchemiops.neighborlist import neighbor_list
import torch
from torch import Tensor


def build_neighbor_list(
    positions: Tensor,
    r_max: float,
    cell: Tensor,
    pbc: Tensor,
    half_fill: bool = True,
    max_neighbors: int | None = None,
    method: str = "cell_list",
):
    """
    Wrap Toolkit-Ops neighbor list to keep a single import site.
    Returns nlist (2, num_pairs), shifts (num_pairs,3).
    """
    nlist, _, shifts = neighbor_list(
        positions,
        float(r_max),
        cell=cell,
        pbc=pbc,
        return_neighbor_list=True,
        half_fill=half_fill,
        max_neighbors=max_neighbors,
        method=method,
    )
    return nlist, shifts
