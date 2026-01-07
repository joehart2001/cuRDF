from collections.abc import Iterable
from typing import Sequence

import numpy as np
from pathlib import Path
from tqdm import tqdm, TqdmWarning
import warnings

warnings.filterwarnings("ignore", category=TqdmWarning)

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
    species_a: str,
    species_b: str | None = None,
    selection: str | None = None,
    selection_b: str | None = None,
    atom_types_map: dict | None = None,
    index=None,
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
    species_a/species_b: required element names for groups A/B. If species_b omitted, uses same-species.
    """
    if mda is None:
        raise ImportError("MDAnalysis must be installed for rdf_from_mdanalysis")
    if torch_dtype is None:
        import torch
        torch_dtype = torch.float32
    warnings.filterwarnings(
        "ignore",
        message="DCDReader currently makes independent timesteps",
        category=DeprecationWarning,
    )

    # Ensure atom names are present
    has_names = hasattr(universe.atoms, "names")
    if not has_names:
        if atom_types_map:
            try:
                types = universe.atoms.types
            except Exception:
                raise ValueError("Topology lacks names and types; cannot map species.")
            mapped = []
            for t in types:
                key_int = None
                try:
                    key_int = int(t)
                except Exception:
                    pass
                name = None
                if key_int is not None and key_int in atom_types_map:
                    name = atom_types_map[key_int]
                elif str(t) in atom_types_map:
                    name = atom_types_map[str(t)]
                if name is None:
                    raise ValueError(f"No mapping in atom_types_map for type '{t}'")
                mapped.append(name)
            universe.add_TopologyAttr("name", mapped)
        elif species_b is None or species_b == species_a:
            universe.add_TopologyAttr("name", [species_a] * len(universe.atoms))
        else:
            raise ValueError("Topology lacks atom names; provide --atom-types mapping for cross-species selections.")

    ag_a = universe.select_atoms(f"name {species_a}")
    ag_b = universe.select_atoms(f"name {species_b}") if species_b is not None else ag_a
    if len(ag_a) == 0:
        raise ValueError(f"No atoms found for species_a='{species_a}' (check atom names or --atom-types mapping)")
    if len(ag_b) == 0:
        raise ValueError(f"No atoms found for species_b='{species_b or species_a}' (check atom names or --atom-types mapping)")
    if wrap_positions and mda_wrap is not None:
        ag_wrap = ag_a if selection_b is None else (ag_a | ag_b)
        universe.trajectory.add_transformations(mda_wrap(ag_wrap, compound="atoms"))

    same_species = len(ag_a) == len(ag_b) and ag_a is ag_b

    def frames():
        traj = universe.trajectory[index] if index is not None else universe.trajectory
        for ts in tqdm(traj, desc="Frames (MDAnalysis)", unit="frame"):
            cell = _mdanalysis_cell_matrix(ts.dimensions)
            if same_species:
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

    if not same_species and half_fill:
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


def rdf(
    obj,
    species_a: str,
    species_b: str | None = None,
    index=None,
    atom_types_map: dict | None = None,
    outdir=None,
    output: str | None = None,
    **kwargs,
):
    """
    Unified entry point:
      - ASE: pass an Atoms or iterable of Atoms
      - MDAnalysis: pass a Universe

    Additional kwargs are forwarded to rdf_from_ase or rdf_from_mdanalysis.
    If output is provided, saves bins/gr to that path (npz). If outdir is provided, saves to outdir/rdf.npz.
    """
    # Lazy imports to avoid hard deps
    if hasattr(obj, "get_positions"):  # ASE Atoms
        bins, gr = rdf_from_ase(
            obj,
            species_a=species_a,
            species_b=species_b,
            index=index,
            atom_types_map=atom_types_map,
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
            **kwargs,
        )
        _maybe_save(outdir, output, bins, gr)
        return bins, gr
    raise TypeError("rdf() expects an ASE Atoms/trajectory or an MDAnalysis Universe")


def _maybe_save(outdir, output, bins, gr):
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


def rdf_from_ase(
    atoms_or_trajectory,
    selection: Sequence[int] | None = None,
    selection_b: Sequence[int] | None = None,
    species_a: str | None = None,
    species_b: str | None = None,
    index=None,
    atom_types_map: dict | None = None,
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
    
    selection/selection_b: index lists for group A and group B (cross-species).
    species_a/species_b: element symbols for group A/B (if provided, override selection indices).
    With only one group provided, computes same-species RDF.
    """
    if Atoms is None:
        raise ImportError("ASE must be installed for rdf_from_ase")
    if torch_dtype is None:
        import torch
        torch_dtype = torch.float32

    def _frames_iter():
        if hasattr(atoms_or_trajectory, "get_positions"):
            iterable = (atoms_or_trajectory,)
        elif index is not None:
            iterable = atoms_or_trajectory[index]
        elif isinstance(atoms_or_trajectory, Iterable):
            iterable = atoms_or_trajectory
        else:
            raise TypeError("atoms_or_trajectory must be ASE Atoms or iterable of Atoms")

        for frame in tqdm(iterable, desc="Frames (ASE)", unit="frame"):
            if not hasattr(frame, "get_positions"):
                raise TypeError("Each frame must be ASE Atoms")
            n_atoms = len(frame)
            symbols = frame.get_chemical_symbols()

            if atom_types_map:
                # map numeric types to element names if provided
                types = frame.get_array("numbers", int, None)
                name_list = []
                for t in types:
                    name = atom_types_map.get(t) or atom_types_map.get(str(t))
                    if name is None:
                        raise ValueError(f"No mapping in atom_types_map for type '{t}'")
                    name_list.append(name)
                symbols = name_list

            if species_a is not None:
                idx_a = np.where(np.array(symbols) == species_a)[0]
            else:
                idx_a = _extract_selection_indices(selection, n_atoms)

            if species_b is not None:
                idx_b = np.where(np.array(symbols) == species_b)[0]
            elif selection_b is not None:
                idx_b = _extract_selection_indices(selection_b, n_atoms)
            else:
                idx_b = idx_a

            if len(idx_a) == 0:
                raise ValueError(f"No atoms found for species/selection A ({species_a or selection})")
            if len(idx_b) == 0:
                raise ValueError(f"No atoms found for species/selection B ({species_b or selection_b or selection})")

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

    same_species = (
        (species_b is None or species_b == species_a)
        and selection_b is None
    )
    if not same_species and half_fill:
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
