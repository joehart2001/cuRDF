from collections.abc import Iterable
from typing import Sequence

import numpy as np
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
    Convert MDAnalysis unit cell representation to a 3x3 matrix.

    Parameters
    ----------
    dimensions
        Array-like with ``[a, b, c, alpha, beta, gamma]`` in MDAnalysis format.

    Returns
    -------
    np.ndarray
        Cell matrix shaped ``(3, 3)``.
    """
    if triclinic_vectors is None:
        raise ImportError("MDAnalysis not available")
    return np.array(triclinic_vectors(dimensions), dtype=np.float32)


def rdf_from_mdanalysis(
    universe,
    species_a: str,
    species_b: str | None = None,
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
    method: str = "cell_list",
):
    """
    Compute g(r) from an MDAnalysis Universe across trajectory frames.

    Parameters
    ----------
    universe
        MDAnalysis ``Universe`` containing topology and trajectory.
    species_a
        Element name for group A.
    species_b
        Optional element name for group B; defaults to group A.
    selection_b
        Optional selection string for group B; defaults to species-based selection.
    atom_types_map
        Optional mapping for numeric atom types to element names.
    index
        Optional trajectory index/selector.
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
    wrap_positions
        Whether to wrap positions into the unit cell before counting.
    method
        Neighbor-list method name (e.g., ``"cell_list"`` or ``"naive"``).

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Bin centers and g(r) arrays on CPU.
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
        method=method,
    )


def _extract_selection_indices(selection: Sequence[int] | None, n_atoms: int):
    """
    Normalize an index selection to a 1D numpy array.

    Parameters
    ----------
    selection
        Optional sequence of indices.
    n_atoms
        Total atom count used for bounds checking.

    Returns
    -------
    np.ndarray
        Validated index array.
    """
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
    method: str = "cell_list",
    wrap_positions: bool = True,
):
    """
    Compute g(r) from an ASE ``Atoms`` or iterable of ``Atoms``.

    Parameters
    ----------
    atoms_or_trajectory
        ASE ``Atoms`` or iterable of ``Atoms`` frames.
    selection
        Optional indices for group A.
    selection_b
        Optional indices for group B; defaults to selection.
    species_a
        Optional element name for group A; overrides selections.
    species_b
        Optional element name for group B; overrides selections.
    index
        Optional trajectory slice/index.
    atom_types_map
        Optional mapping for numeric atom types to element names.
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
    wrap_positions
        Whether to wrap positions into the unit cell before counting.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Bin centers and g(r) arrays on CPU.
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
        method=method,
    )
