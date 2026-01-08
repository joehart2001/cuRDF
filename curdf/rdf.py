import math
from collections.abc import Iterable as CollIterable
from pathlib import Path
from typing import Iterable

import numpy as np
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

    Parameters
    ----------
    counts
        Histogram tensor updated in place.
    positions
        Cartesian coordinates shaped ``(N, 3)``.
    cell
        Cell tensor shaped ``(1, 3, 3)``.
    pbc
        Periodic boundary flags tensor shaped ``(1, 3)``.
    edges
        Histogram bin edges tensor shaped ``(nbins + 1,)``.
    r_min
        Minimum distance included in the histogram.
    r_max
        Maximum distance included in the histogram.
    half_fill
        Whether to build unique pairs (same-species mode).
    max_neighbors
        Optional neighbor-list cap forwarded to Toolkit-Ops.
    method
        Neighbor-list method name passed to Toolkit-Ops.
    group_a_mask
        Optional boolean mask selecting source atoms.
    group_b_mask
        Optional boolean mask selecting target atoms.

    Returns
    -------
    float
        Normalization factor ``n_group_a * rho_group_b`` so callers can combine frames.
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
    """
    Convert accumulated counts into g(r).

    Parameters
    ----------
    counts
        Histogram counts tensor.
    edges
        Bin edges tensor shaped ``(nbins + 1,)``.
    total_norm
        Accumulated normalization factor from frames.
    half_fill
        Whether pairs were unique (affects pair factor).
    cross_mode
        Whether cross-species mode was used (affects pair factor).

    Returns
    -------
    tuple[Tensor, Tensor]
        Bin centers tensor and g(r) tensor.
    """
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

    Parameters
    ----------
    positions
        Array-like Cartesian coordinates shaped ``(N, 3)``.
    cell
        Cell matrix shaped ``(3, 3)`` (triclinic allowed).
    pbc
        Iterable of three booleans for periodicity along each axis.
    r_min
        Minimum distance included in the histogram.
    r_max
        Maximum distance included in the histogram.
    nbins
        Number of histogram bins.
    device
        Torch device string or object used for computation.
    torch_dtype
        Torch dtype used for tensors.
    half_fill
        ``True`` for identical-species mode (unique pairs), ``False`` for ordered pairs.
    max_neighbors
        Optional neighbor-list cap forwarded to Toolkit-Ops.
    method
        Neighbor-list method name (e.g., ``"cell_list"`` or ``"naive"``).
    group_a_indices
        Optional indices for group A (cross-species mode).
    group_b_indices
        Optional indices for group B (cross-species mode); defaults to group A when omitted.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Bin centers and g(r) arrays on CPU.
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
    method: str = "cell_list",
):
    """
    Accumulate g(r) over multiple frames.

    Parameters
    ----------
    frames
        Iterable yielding dicts with keys ``positions``, ``cell``, ``pbc``, and optional ``group_a_mask``/``group_b_mask``.
    r_min
        Minimum distance included in the histogram.
    r_max
        Maximum distance included in the histogram.
    nbins
        Number of histogram bins.
    device
        Torch device string or object used for computation.
    torch_dtype
        Torch dtype used for tensors.
    half_fill
        ``True`` for identical-species mode (unique pairs), ``False`` for ordered pairs.
    max_neighbors
        Optional neighbor-list cap forwarded to Toolkit-Ops.
    method
        Neighbor-list method name (e.g., ``"cell_list"`` or ``"naive"``).

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Bin centers and g(r) arrays on CPU.
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


def rdf(
    obj,
    species_a: str,
    species_b: str | None = None,
    index=None,
    atom_types_map: dict | None = None,
    method: str = "cell_list",
    outdir=None,
    output: str | None = None,
    **kwargs,
):
    """
    Unified user entry point for g(r) computation.

    Parameters
    ----------
    obj
        ASE ``Atoms``, iterable of ``Atoms`` (e.g., a list trajectory), or MDAnalysis ``Universe``.
    species_a
        Element name for group A.
    species_b
        Optional element name for group B (defaults to group A).
    index
        Optional index/selector forwarded to source-specific readers.
    atom_types_map
        Optional mapping for numeric atom types to element names.
    method
        Neighbor-list method name (e.g., ``"cell_list"`` or ``"naive"``).
    outdir
        Optional directory to save ``rdf.npz``.
    output
        Optional explicit output path (``.npz``, ``.csv``, ``.json``, ``.pkl``).
    **kwargs
        Additional parameters forwarded to source-specific readers.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Bin centers and g(r) arrays on CPU.
    """
    # Lazy imports to avoid hard deps
    from .adapters import rdf_from_ase, rdf_from_mdanalysis

    if hasattr(obj, "get_positions"):  # ASE Atoms
        bins, gr = rdf_from_ase(
            obj,
            species_a=species_a,
            species_b=species_b,
            index=index,
            atom_types_map=atom_types_map,
            method=method,
            **kwargs,
        )
        _maybe_save(outdir, output, bins, gr)
        return bins, gr
    if isinstance(obj, CollIterable) and all(hasattr(f, "get_positions") for f in obj):
        bins, gr = rdf_from_ase(
            obj,
            species_a=species_a,
            species_b=species_b,
            index=index,
            atom_types_map=atom_types_map,
            method=method,
            **kwargs,
        )
        _maybe_save(outdir, output, bins, gr)
        return bins, gr
    # MDAnalysis Universe duck check
    if hasattr(obj, "trajectory") and hasattr(obj, "select_atoms"):
        bins, gr = rdf_from_mdanalysis(
            obj,
            species_a=species_a,
            species_b=species_b,
            index=index,
            atom_types_map=atom_types_map,
            method=method,
            **kwargs,
        )
        _maybe_save(outdir, output, bins, gr)
        return bins, gr
    raise TypeError("rdf() expects an ASE Atoms/trajectory or an MDAnalysis Universe")


def _maybe_save(outdir, output, bins, gr):
    """
    Persist RDF results to disk when requested.

    Parameters
    ----------
    outdir
        Optional output directory; saves ``rdf.npz`` when provided.
    output
        Optional explicit output path.
    bins
        Bin centers array.
    gr
        g(r) array matching ``bins``.
    """
    target = None
    if output:
        target = Path(output)
        target.parent.mkdir(parents=True, exist_ok=True)
    elif outdir:
        path = Path(outdir)
        path.mkdir(parents=True, exist_ok=True)
        target = path / "rdf.npz"
    else:
        return

    suffix = target.suffix.lower()
    if suffix == ".npz":
        np.savez(target, bins=bins, gr=gr)
    elif suffix in {".csv", ".tsv"}:
        sep = "," if suffix == ".csv" else "\t"
        import pandas as pd

        df = pd.DataFrame({"r": bins, "g_r": gr})
        df.to_csv(target, sep=sep, index=False)
    elif suffix in {".json"}:
        import json

        with open(target, "w") as f:
            json.dump({"bins": bins.tolist(), "g_r": gr.tolist()}, f)
    elif suffix in {".pkl", ".pickle"}:
        import pickle

        with open(target, "wb") as f:
            pickle.dump({"bins": bins, "g_r": gr}, f, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        np.savez(target, bins=bins, gr=gr)
