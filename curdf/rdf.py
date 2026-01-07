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
    max_neighbors: int | None,
    method: str,
    group_a_mask: Tensor | None = None,
    group_b_mask: Tensor | None = None,
) -> float:
    """
    Accumulate pair counts for one frame.
    Returns normalization factor (n_group_a * rho_group_b) so the caller can normalize after multiple frames.
    """
    dr = (r_max - r_min) / (len(edges) - 1)

    nlist, shifts = build_neighbor_list(
        positions,
        r_max,
        cell=cell,
        pbc=pbc,
        half_fill=half_fill,
        max_neighbors=max_neighbors,
        method=method,
    )

    src = nlist[0].to(torch.int64)
    tgt = nlist[1].to(torch.int64)

    shift_cart = (shifts.to(positions.dtype) @ cell[0])
    dr_vec = (positions[tgt] + shift_cart) - positions[src]
    dist = torch.linalg.norm(dr_vec, dim=1)

    valid = (dist >= r_min) & (dist < r_max)
    if group_a_mask is not None and group_b_mask is not None:
        src_mask = group_a_mask[src]
        tgt_mask = group_b_mask[tgt]
        valid = valid & src_mask & tgt_mask
    dist = dist[valid]

    bin_idx = torch.floor((dist - r_min) / dr).to(torch.int64)
    bin_idx = torch.clamp(bin_idx, 0, counts.numel() - 1)
    counts.scatter_add_(0, bin_idx, torch.ones_like(bin_idx, dtype=torch.int64))

    volume = cell_volume(cell)
    if group_a_mask is not None and group_b_mask is not None:
        n_a = group_a_mask.sum().item()
        n_b = group_b_mask.sum().item()
        norm_factor = n_a * (n_b / volume)
    else:
        n_atoms = positions.shape[0]
        norm_factor = n_atoms * (n_atoms / volume)  # n_atoms * rho
    return norm_factor


def _finalize_gr(
    counts: Tensor,
    edges: Tensor,
    total_norm: float,
    half_fill: bool,
    cross_mode: bool,
) -> tuple[Tensor, Tensor]:
    r1 = edges[:-1]
    r2 = edges[1:]
    shell_vol = (4.0 / 3.0) * math.pi * (r2**3 - r1**3)
    pair_factor = 1.0 if cross_mode else (2.0 if half_fill else 1.0)

    if total_norm == 0:
        raise ValueError("Total normalization is zero; no frames processed or selections/species empty?")
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
    max_neighbors: int | None = None,
    method: str = "cell_list",
    group_a_indices=None,
    group_b_indices=None,
):
    """
    Compute g(r) for a single frame of positions.

    Args:
        positions: array-like (N,3)
        cell: (3,3) cell matrix (triclinic allowed)
        pbc: iterable of 3 booleans
        r_min/r_max/nbins: histogram parameters
        half_fill: True for identical species (unique pairs); False for ordered pairs
        max_neighbors: passed to Toolkit-Ops neighbor list (None for library default)
        method: neighbor list method ("cell_list" or "naive")
        group_a_indices/group_b_indices: optional index lists for cross-species RDF.
            If provided, counts pairs with src in A and tgt in B. When both are None,
            uses all atoms (identical-species mode).
    """
    device = torch.device(device)
    pos_t = torch.as_tensor(positions, device=device, dtype=torch_dtype)
    if pos_t.ndim != 2 or pos_t.shape[1] != 3:
        raise ValueError(f"positions must be (N,3); got {tuple(pos_t.shape)}")

    cell_t = cell_tensor(cell, device=device, dtype=torch_dtype)
    pbc_t = pbc_tensor(pbc, device=device)

    edges = torch.linspace(r_min, r_max, nbins + 1, device=device, dtype=torch_dtype)
    counts = torch.zeros(nbins, device=device, dtype=torch.int64)

    group_a_mask = group_b_mask = None
    if group_a_indices is not None:
        group_a_mask = torch.zeros(pos_t.shape[0], device=device, dtype=torch.bool)
        group_a_mask[torch.as_tensor(group_a_indices, device=device, dtype=torch.long)] = True
    if group_b_indices is not None:
        group_b_mask = torch.zeros(pos_t.shape[0], device=device, dtype=torch.bool)
        group_b_mask[torch.as_tensor(group_b_indices, device=device, dtype=torch.long)] = True
    elif group_a_mask is not None:
        group_b_mask = group_a_mask

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
        method=method,
        group_a_mask=group_a_mask,
        group_b_mask=group_b_mask,
    )

    cross_mode = group_a_mask is not None and group_b_mask is not None and not torch.equal(
        group_a_mask, group_b_mask
    )
    centers, g_r = _finalize_gr(
        counts, edges, total_norm, half_fill=half_fill, cross_mode=cross_mode
    )
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
    frames: iterable yielding dicts with keys positions, cell, pbc, and optional group_a_mask/group_b_mask
    """
    device = torch.device(device)
    edges = torch.linspace(r_min, r_max, nbins + 1, device=device, dtype=torch_dtype)
    counts = torch.zeros(nbins, device=device, dtype=torch.int64)
    total_norm = 0.0

    cross_flag = False

    for frame in frames:
        pos_t = torch.as_tensor(frame["positions"], device=device, dtype=torch_dtype)
        cell_t = cell_tensor(frame["cell"], device=device, dtype=torch_dtype)
        pbc_t = pbc_tensor(frame["pbc"], device=device)
        group_a_mask = frame.get("group_a_mask")
        group_b_mask = frame.get("group_b_mask")
        if group_a_mask is not None:
            group_a_mask = torch.as_tensor(group_a_mask, device=device, dtype=torch.bool)
        if group_b_mask is not None:
            group_b_mask = torch.as_tensor(group_b_mask, device=device, dtype=torch.bool)
        elif group_a_mask is not None:
            group_b_mask = group_a_mask
        if group_a_mask is not None and group_b_mask is not None and not torch.equal(
            group_a_mask, group_b_mask
        ):
            cross_flag = True

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
        method=method,
        group_a_mask=group_a_mask,
        group_b_mask=group_b_mask,
    )
        total_norm += norm

    centers, g_r = _finalize_gr(
        counts,
        edges,
        total_norm,
        half_fill=half_fill,
        cross_mode=cross_flag,
    )
    return centers.cpu().numpy(), g_r.cpu().numpy()
