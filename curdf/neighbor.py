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
    Build a neighbor list via Toolkit-Ops.

    Parameters
    ----------
    positions
        Cartesian coordinates shaped ``(N, 3)``.
    r_max
        Cutoff radius for neighbor search.
    cell
        Cell tensor shaped ``(1, 3, 3)``.
    pbc
        Boolean PBC tensor shaped ``(1, 3)``.
    half_fill
        ``True`` for unique pairs, ``False`` for ordered pairs.
    max_neighbors
        Optional neighbor cap forwarded to Toolkit-Ops.
    method
        Neighbor-list method name (e.g., ``"cell_list"`` or ``"naive"``).

    Returns
    -------
    tuple[Tensor, Tensor]
        Neighbor index tensor shaped ``(2, num_pairs)`` and shift vectors shaped ``(num_pairs, 3)``.
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
