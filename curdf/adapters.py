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
    """
    if mda is None:
        raise ImportError("MDAnalysis must be installed for rdf_from_mdanalysis")
    if torch_dtype is None:
        import torch
        torch_dtype = torch.float32

    ag = universe.select_atoms(selection)
    if wrap_positions and mda_wrap is not None:
        universe.trajectory.add_transformations(mda_wrap(ag, compound="atoms"))

    def frames():
        for ts in universe.trajectory:
            yield {
                "positions": ag.positions.astype(np.float32, copy=False),
                "cell": _mdanalysis_cell_matrix(ts.dimensions),
                "pbc": (True, True, True),
            }

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
    selection: list/array of atom indices to include.
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
            idx = _extract_selection_indices(selection, n_atoms)
            pos = frame.get_positions(wrap=wrap_positions)[idx]
            cell = np.array(frame.get_cell().array, dtype=np.float32)
            pbc = tuple(bool(x) for x in frame.get_pbc())

            yield {"positions": pos.astype(np.float32, copy=False), "cell": cell, "pbc": pbc}

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
