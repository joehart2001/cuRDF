from collections.abc import Iterable
from typing import Sequence

import numpy as np

try:
    import MDAnalysis as mda
    from MDAnalysis.lib.mdamath import triclinic_vectors
    from MDAnalysis.transformations import wrap as mda_wrap
except ImportError:
    mda = None
    triclinic_vectors = None
    mda_wrap = None

try:
    from ase import Atoms
except ImportError:
    Atoms = None

from .rdf import accumulate_rdf


def _mdanalysis_cell_matrix(dimensions):
    """
    MDAnalysis gives [a, b, c, alpha, beta, gamma]; convert to 3x3.
    """
    if triclinic_vectors is None:
        raise ImportError("MDAnalysis not available")
    return np.array(triclinic_vectors(dimensions), dtype=np.float32)


def rdf_from_mdanalysis(
    universe,
    selection: str = "all",
    selection_b: str | None = None,
    r_min: float = 1.0,
    r_max: float = 6.0,
    nbins: int = 100,
    device="cuda",
    torch_dtype=None,
    half_fill: bool = True,
    max_neighbors: int = 2048,
    wrap_positions: bool = True,
):
    """
    Compute g(r) from an MDAnalysis Universe across all trajectory frames.
    selection_b: optional second selection for cross-species RDF (A in selection, B in selection_b).
    """
    if mda is None:
        raise ImportError("MDAnalysis must be installed for rdf_from_mdanalysis")
    if torch_dtype is None:
        import torch
        torch_dtype = torch.float32

    ag_a = universe.select_atoms(selection)
    ag_b = universe.select_atoms(selection_b) if selection_b is not None else ag_a
    if wrap_positions and mda_wrap is not None:
        ag_wrap = ag_a if selection_b is None else (ag_a | ag_b)
        universe.trajectory.add_transformations(mda_wrap(ag_wrap, compound="atoms"))

    def frames():
        for ts in universe.trajectory:
            cell = _mdanalysis_cell_matrix(ts.dimensions)
            if selection_b is None:
                yield {
                    "positions": ag_a.positions.astype(np.float32, copy=False),
                    "cell": cell,
                    "pbc": (True, True, True),
                }
            else:
                pos_a = ag_a.positions.astype(np.float32, copy=False)
                pos_b = ag_b.positions.astype(np.float32, copy=False)
                pos = np.concatenate([pos_a, pos_b], axis=0)
                group_a_mask = np.zeros(len(pos), dtype=bool)
                group_b_mask = np.zeros(len(pos), dtype=bool)
                group_a_mask[: len(pos_a)] = True
                group_b_mask[len(pos_a) :] = True
                yield {
                    "positions": pos,
                    "cell": cell,
                    "pbc": (True, True, True),
                    "group_a_mask": group_a_mask,
                    "group_b_mask": group_b_mask,
                }

    if selection_b is not None and half_fill:
        half_fill = False  # cross-species -> ordered pairs

    return accumulate_rdf(
        frames(),
        r_min=r_min,
        r_max=r_max,
        nbins=nbins,
        device=device,
        torch_dtype=torch_dtype,
        half_fill=half_fill,
        max_neighbors=max_neighbors,
    )


def _extract_selection_indices(selection: Sequence[int] | None, n_atoms: int):
    if selection is None:
        return np.arange(n_atoms)
    idx = np.asarray(selection, dtype=int)
    if idx.ndim != 1:
        raise ValueError("selection indices must be 1D")
    if idx.min(initial=0) < 0 or idx.max(initial=0) >= n_atoms:
        raise ValueError("selection indices out of bounds")
    return idx


def rdf_from_ase(
    atoms_or_trajectory,
    selection: Sequence[int] | None = None,
    selection_b: Sequence[int] | None = None,
    r_min: float = 1.0,
    r_max: float = 6.0,
    nbins: int = 100,
    device="cuda",
    torch_dtype=None,
    half_fill: bool = True,
    max_neighbors: int = 2048,
    wrap_positions: bool = True,
):
    """
    Compute g(r) from an ASE Atoms or iterable of Atoms (trajectory).
    selection/selection_b: index lists for group A and group B (cross-species). With only selection provided, computes A–A.
    """
    if Atoms is None:
        raise ImportError("ASE must be installed for rdf_from_ase")
    if torch_dtype is None:
        import torch
        torch_dtype = torch.float32

    def _frames_iter():
        if hasattr(atoms_or_trajectory, "get_positions"):
            iterable = (atoms_or_trajectory,)
        elif isinstance(atoms_or_trajectory, Iterable):
            iterable = atoms_or_trajectory
        else:
            raise TypeError("atoms_or_trajectory must be ASE Atoms or iterable of Atoms")

        for frame in iterable:
            if not hasattr(frame, "get_positions"):
                raise TypeError("Each frame must be ASE Atoms")
            n_atoms = len(frame)
            idx_a = _extract_selection_indices(selection, n_atoms)
            idx_b = _extract_selection_indices(selection_b, n_atoms) if selection_b is not None else idx_a
            pos_all = frame.get_positions(wrap=wrap_positions)
            pos_a = pos_all[idx_a]
            pos_b = pos_all[idx_b]
            pos = np.concatenate([pos_a, pos_b], axis=0)
            cell = np.array(frame.get_cell().array, dtype=np.float32)
            pbc = tuple(bool(x) for x in frame.get_pbc())
            group_a_mask = np.zeros(len(pos), dtype=bool)
            group_b_mask = np.zeros(len(pos), dtype=bool)
            group_a_mask[: len(pos_a)] = True
            group_b_mask[len(pos_a) :] = True

            yield {
                "positions": pos.astype(np.float32, copy=False),
                "cell": cell,
                "pbc": pbc,
                "group_a_mask": group_a_mask,
                "group_b_mask": group_b_mask,
            }

    if selection_b is not None and half_fill:
        half_fill = False  # cross-species -> ordered pairs

    return accumulate_rdf(
        _frames_iter(),
        r_min=r_min,
        r_max=r_max,
        nbins=nbins,
        device=device,
        torch_dtype=torch_dtype,
        half_fill=half_fill,
        max_neighbors=max_neighbors,
    )
