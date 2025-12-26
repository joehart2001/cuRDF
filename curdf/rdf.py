import math
from typing import Iterable

import torch
from torch import Tensor

from .cell import cell_tensor, cell_volume, pbc_tensor
from .neighbor import build_neighbor_list


def _update_counts(
    counts: Tensor,
    positions: Tensor,
    cell: Tensor,
    pbc: Tensor,
    edges: Tensor,
    r_min: float,
    r_max: float,
    half_fill: bool,
    max_neighbors: int,
) -> float:
    """
    Accumulate pair counts for one frame.
    Returns rho for this frame so the caller can normalize after multiple frames.
    """
    dr = (r_max - r_min) / (len(edges) - 1)

    nlist, shifts = build_neighbor_list(
        positions,
        r_max,
        cell=cell,
        pbc=pbc,
        half_fill=half_fill,
        max_neighbors=max_neighbors,
    )

    src = nlist[0].to(torch.int64)
    tgt = nlist[1].to(torch.int64)

    shift_cart = (shifts.to(positions.dtype) @ cell[0])
    dr_vec = (positions[tgt] + shift_cart) - positions[src]
    dist = torch.linalg.norm(dr_vec, dim=1)

    valid = (dist >= r_min) & (dist < r_max)
    dist = dist[valid]

    bin_idx = torch.floor((dist - r_min) / dr).to(torch.int64)
    bin_idx = torch.clamp(bin_idx, 0, counts.numel() - 1)
    counts.scatter_add_(0, bin_idx, torch.ones_like(bin_idx, dtype=torch.int64))

    volume = cell_volume(cell)
    n_atoms = positions.shape[0]
    norm_factor = n_atoms * (n_atoms / volume)  # n_atoms * rho
    return norm_factor


def _finalize_gr(counts: Tensor, edges: Tensor, total_norm: float, half_fill: bool) -> tuple[Tensor, Tensor]:
    r1 = edges[:-1]
    r2 = edges[1:]
    shell_vol = (4.0 / 3.0) * math.pi * (r2**3 - r1**3)
    pair_factor = 2.0 if half_fill else 1.0

    if total_norm == 0:
        raise ValueError("Total normalization is zero; no frames processed?")
    g_r = (pair_factor * counts.to(r1.dtype)) / (shell_vol * total_norm)
    centers = (edges[:-1] + edges[1:]) * 0.5
    return centers, g_r


@torch.no_grad()
def compute_rdf(
    positions,
    cell,
    pbc=(True, True, True),
    r_min: float = 1.0,
    r_max: float = 6.0,
    nbins: int = 100,
    device: str | torch.device = "cuda",
    torch_dtype: torch.dtype = torch.float32,
    half_fill: bool = True,
    max_neighbors: int = 2048,
):
    """
    Compute g(r) for a single frame of positions.

    Args:
        positions: array-like (N,3)
        cell: (3,3) cell matrix (triclinic allowed)
        pbc: iterable of 3 booleans
        r_min/r_max/nbins: histogram parameters
        half_fill: True for identical species (unique pairs); False for ordered pairs
        max_neighbors: passed to Toolkit-Ops neighbor list
    """
    device = torch.device(device)
    pos_t = torch.as_tensor(positions, device=device, dtype=torch_dtype)
    if pos_t.ndim != 2 or pos_t.shape[1] != 3:
        raise ValueError(f"positions must be (N,3); got {tuple(pos_t.shape)}")

    cell_t = cell_tensor(cell, device=device, dtype=torch_dtype)
    pbc_t = pbc_tensor(pbc, device=device)

    edges = torch.linspace(r_min, r_max, nbins + 1, device=device, dtype=torch_dtype)
    counts = torch.zeros(nbins, device=device, dtype=torch.int64)

    total_norm = _update_counts(
        counts,
        pos_t,
        cell=cell_t,
        pbc=pbc_t,
        edges=edges,
        r_min=r_min,
        r_max=r_max,
        half_fill=half_fill,
        max_neighbors=max_neighbors,
    )

    centers, g_r = _finalize_gr(counts, edges, total_norm, half_fill=half_fill)
    return centers.cpu().numpy(), g_r.cpu().numpy()


@torch.no_grad()
def accumulate_rdf(
    frames: Iterable[dict],
    r_min: float,
    r_max: float,
    nbins: int,
    device: str | torch.device,
    torch_dtype: torch.dtype,
    half_fill: bool,
    max_neighbors: int,
):
    """
    General accumulator for multiple frames.
    frames: iterable yielding dicts with keys positions, cell, pbc
    """
    device = torch.device(device)
    edges = torch.linspace(r_min, r_max, nbins + 1, device=device, dtype=torch_dtype)
    counts = torch.zeros(nbins, device=device, dtype=torch.int64)
    total_norm = 0.0

    for frame in frames:
        pos_t = torch.as_tensor(frame["positions"], device=device, dtype=torch_dtype)
        cell_t = cell_tensor(frame["cell"], device=device, dtype=torch_dtype)
        pbc_t = pbc_tensor(frame["pbc"], device=device)

        norm = _update_counts(
            counts,
            pos_t,
            cell=cell_t,
            pbc=pbc_t,
            edges=edges,
            r_min=r_min,
            r_max=r_max,
            half_fill=half_fill,
            max_neighbors=max_neighbors,
        )
        total_norm += norm

    centers, g_r = _finalize_gr(counts, edges, total_norm, half_fill=half_fill)
    return centers.cpu().numpy(), g_r.cpu().numpy()
