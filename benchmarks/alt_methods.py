import time
from typing import Sequence

import numpy as np


def _bin_distances(dist: np.ndarray, edges: np.ndarray) -> None:
    """Accumulate histogram counts for timing parity with RDF calculations."""
    _ = np.histogram(dist, bins=edges)


def bench_mdanalysis(frames: Sequence[np.ndarray], cell: np.ndarray, r_min: float, r_max: float, nbins: int):
    """CPU reference using MDAnalysis; returns elapsed seconds or None if unavailable."""
    try:
        import MDAnalysis as mda
        from MDAnalysis.analysis.rdf import InterRDF
        from MDAnalysis.coordinates.memory import MemoryReader
        from MDAnalysis.transformations import wrap
    except ImportError:
        print("Skipping MDAnalysis: package not installed")
        return None

    coords = np.stack(frames)
    dims = np.tile(
        np.array([cell[0, 0], cell[1, 1], cell[2, 2], 90.0, 90.0, 90.0], dtype=np.float32),
        (coords.shape[0], 1),
    )

    start = time.perf_counter()
    u = mda.Universe.empty(coords.shape[1], trajectory=True)
    u.add_TopologyAttr("name", ["C"] * coords.shape[1])
    u.trajectory = MemoryReader(coords, dimensions=dims)
    atoms = u.select_atoms("all")
    u.trajectory.add_transformations(wrap(atoms, compound="atoms"))

    rdf = InterRDF(atoms, atoms, range=(r_min, r_max), nbins=nbins)
    rdf.run()
    return time.perf_counter() - start


def bench_matscipy(frames: Sequence[np.ndarray], cell: np.ndarray, edges: np.ndarray, r_min: float, r_max: float):
    """Bench matscipy.neighbours.neighbour_list; returns elapsed seconds or None if unavailable."""
    try:
        from matscipy.neighbours import neighbour_list
    except ImportError:
        print("Skipping matscipy: package not installed")
        return None

    pbc = np.array([True, True, True], dtype=bool)
    start = time.perf_counter()
    for pos in frames:
        i, j, S = neighbour_list("ijS", positions=pos, cell=cell, pbc=pbc, cutoff=r_max)
        mask = i < j  # avoid double counting
        i = i[mask]
        j = j[mask]
        S = S[mask]
        dr = pos[j] + S @ cell - pos[i]
        dist = np.linalg.norm(dr, axis=1)
        dist = dist[(dist >= r_min) & (dist < r_max)]
        _bin_distances(dist, edges)
    return time.perf_counter() - start


def bench_ase(frames: Sequence[np.ndarray], cell: np.ndarray, edges: np.ndarray, r_min: float, r_max: float):
    """Bench ASE NeighborList; returns elapsed seconds or None if unavailable."""
    try:
        from ase.neighborlist import neighbor_list as ase_neighbor_list
        from ase import Atoms as ASEAtoms
    except ImportError:
        print("Skipping ASE neighbor list: package not installed")
        return None

    start = time.perf_counter()
    for pos in frames:
        atoms = ASEAtoms(positions=pos, cell=cell, pbc=True, symbols=["X"] * pos.shape[0])
        i, j, S = ase_neighbor_list("ijS", atoms, r_max, self_interaction=False)
        mask = i < j
        i = i[mask]
        j = j[mask]
        S = S[mask]
        dr = pos[j] + S @ cell - pos[i]
        dist = np.linalg.norm(dr, axis=1)
        dist = dist[(dist >= r_min) & (dist < r_max)]
        _bin_distances(dist, edges)
    return time.perf_counter() - start


def bench_pymatgen(frames: Sequence[np.ndarray], cell: np.ndarray, edges: np.ndarray, r_min: float, r_max: float):
    """Bench pymatgen neighbor list; returns elapsed seconds or None if unavailable."""
    try:
        from pymatgen.core import Lattice, Structure
    except ImportError:
        print("Skipping pymatgen: package not installed")
        return None

    a, b, c = cell[0, 0], cell[1, 1], cell[2, 2]
    lattice = Lattice.from_parameters(a=a, b=b, c=c, alpha=90, beta=90, gamma=90)
    inv_cell = np.linalg.inv(cell)

    start = time.perf_counter()
    for pos in frames:
        frac = pos @ inv_cell  # fractional coords for pymatgen
        structure = Structure(lattice, ["X"] * pos.shape[0], frac, to_unit_cell=True, coords_are_cartesian=False)
        i, j, images, dist = structure.get_neighbor_list(r=r_max)
        mask = i < j
        dist = dist[mask]
        dist = dist[(dist >= r_min) & (dist < r_max)]
        _bin_distances(dist, edges)
    return time.perf_counter() - start
