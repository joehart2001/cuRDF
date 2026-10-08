# API

## `rdf`: structures and trajectories

```python
from curdf import rdf

bins, gr = rdf(
    obj,
    species_a,
    species_b=None,
    index=None,
    atom_types_map=None,
    pbc=None,
    method="cell_list",
    outdir=None,
    output=None,
    **kwargs,
)
```

`obj` must be an ASE `Atoms`, a materialized list of ASE frames, or an
MDAnalysis `Universe`. Read structure files with ASE or MDAnalysis first.
Use a list rather than a one-shot generator, because the entry point inspects
the iterable before processing it.

The return value is `(bins, gr)`, two one-dimensional NumPy arrays of length
`nbins`. `bins` contains bin centres, not edges. Distances are in angstroms
when inputs use angstroms, and `gr` is dimensionless. No uncertainty array is
returned.

### Entry-point parameters

| Parameter | Default | Meaning |
|---|---|---|
| `species_a` | Required | Group A. ASE uses chemical symbols, MDAnalysis uses atom names. |
| `species_b` | `None` | Omit for a same-species RDF. For cross-species calculations, supply a different species. |
| `index` | `None` | Source-specific frame selection. Prefer `slice(start, stop, step)` for a list or an MDAnalysis trajectory. ASE reader strings such as `"1000::10"` belong in `ase.io.read()`. |
| `atom_types_map` | `None` | Integer or string type keys mapped to names. In MDAnalysis, used only when atom names are missing. In ASE, maps values in the `numbers` array. |
| `pbc` | `None` | Defaults to `(True, True, True)`, independently of the source's flags. Pass three booleans to override. |
| `method` | `"cell_list"` | Toolkit-Ops neighbour-list method. `"naive"` is also supported. |
| `outdir` | `None` | Creates the directory and saves `rdf.npz`. |
| `output` | `None` | Explicit output filename. Takes precedence over `outdir`. |

For same-species MDAnalysis calculations, omit `species_b` so both sides use
the same atom group. For custom groups, use identical or disjoint selections.
Overlapping, unequal groups do not have an overlap correction in the
normalization.

### Calculation parameters

These keyword arguments are forwarded to the ASE or MDAnalysis adapter:

| Parameter | Default | Meaning |
|---|---|---|
| `r_min` | `1.0` | Requested lower histogram edge. Set `0.0` to include distances below 1 angstrom. |
| `r_max` | `6.0` | Upper histogram edge and neighbour cutoff. Distances equal to this edge are excluded. |
| `nbins` | `100` | Number of equal-width radial bins. Must be a positive integer. |
| `r_min_floor` | `1e-6` | Effective lower edge is `max(r_min, r_min_floor)`. |
| `device` | `"cuda"` | Torch device. Production calculations use CUDA and do not automatically fall back to CPU. |
| `torch_dtype` | `None` | Resolves to `torch.float32`. Use `torch.float64` for the tensor calculations if needed. |
| `half_fill` | `True` | Unique pairs for same-species calculations. Cross-species calculations force ordered pairs. |
| `max_neighbors` | `2048` | Neighbour capacity forwarded to Toolkit-Ops. Increase when the backend reports insufficient capacity. |
| `wrap_positions` | `True` | Wrap coordinates when all three axes are periodic. |

Keep all lengths in the same units, provide a nonsingular cell, and require
`r_max > max(r_min, r_min_floor)`. See [Methods](methods.md) for normalization
and boundary conventions.

The ASE adapter reads positions and cells into float32 before tensor conversion.
Selecting `torch.float64` does not restore precision lost in that conversion.
The MDAnalysis adapter likewise converts positions and the cell to float32.

### Frame selection and source behavior

For an ASE list, either select frames before calling `rdf()` or pass a Python
slice through `index`. A single ASE `Atoms` is analysed as one frame. To choose
one frame from a list, use `rdf(frames[k], ...)`.

For MDAnalysis, use `start=None`, `stop=None`, and `step=None` to process the
whole trajectory, or set them to choose a frame range. An explicit `index`
takes precedence over these arguments. Use `index=slice(...)` to retain a
trajectory iterator. `stop` is exclusive.

With wrapping enabled, the MDAnalysis adapter adds a wrapping transformation
to the universe for fully periodic cells. If the universe already has
transformations, use `wrap_positions=False` and manage wrapping in your own
pipeline.

### Output formats

| Extension | Contents |
|---|---|
| `.npz` | NumPy arrays named `bins` and `gr`. |
| `.csv`, `.tsv` | Table columns `r` and `g_r`, without an index. Uses pandas. |
| `.json` | Lists under keys `bins` and `g_r`. |
| `.pkl`, `.pickle` | Dictionary with NumPy arrays under keys `bins` and `g_r`. |

Use an explicit supported extension. An unrecognized extension takes the NumPy
NPZ path and may have `.npz` appended to its name. Files contain the curve,
not a record of the cell, selected frames, or calculation settings. Record those
settings separately for reproducibility.

## `compute_rdf`: one frame of coordinates

```python
import torch
from curdf import compute_rdf

bins, gr = compute_rdf(
    positions,
    cell,
    pbc=(True, True, True),
    r_min=1.0,
    r_max=6.0,
    nbins=100,
    r_min_floor=1e-6,
    device="cuda",
    torch_dtype=torch.float32,
    half_fill=True,
    max_neighbors=None,
    method="cell_list",
    group_a_indices=None,
    group_b_indices=None,
)
```

This lower-level function accepts coordinates shaped `(N, 3)` and a cell matrix
shaped `(3, 3)`, with lattice vectors as rows. It also accepts Torch tensors.
Supply three periodicity flags, not a scalar boolean. It returns the same
`(bins, gr)` arrays as `rdf()` and does not write output files.

Without group indices, all atoms form one group. With `group_a_indices` alone,
group B defaults to A. For a cross-group RDF, supply disjoint index lists for
both groups. An example using groups defined by ASE symbols is:

```python
import numpy as np
from ase.io import read
from curdf import compute_rdf

atoms = read("structure.extxyz")
symbols = np.asarray(atoms.get_chemical_symbols())
bins, gr = compute_rdf(
    atoms.positions,
    atoms.cell.array,
    pbc=tuple(atoms.pbc),
    group_a_indices=np.flatnonzero(symbols == "C"),
    group_b_indices=np.flatnonzero(symbols == "O"),
    r_min=0.0,
    r_max=6.0,
    nbins=200,
)
```

Unlike the adapters, this function defaults `max_neighbors` to `None` and
converts inputs directly to the requested tensor dtype.

## `accumulate_rdf`: a frame iterator

```python
import torch
from curdf.rdf import accumulate_rdf

bins, gr = accumulate_rdf(
    frames,
    r_min=0.0,
    r_max=6.0,
    nbins=200,
    r_min_floor=1e-6,
    device="cuda",
    torch_dtype=torch.float32,
    half_fill=True,
    max_neighbors=2048,
    method="cell_list",
)
```

This building block is available from `curdf.rdf`, rather than the top-level
package. Each frame is a dictionary with `positions`, `cell`, and `pbc`.
Optional Boolean arrays `group_a_mask` and `group_b_mask` select groups. If
only A's mask is supplied, B defaults to A. Use the same kind of grouping across
all frames, and identical or disjoint groups within each frame.

All arguments shown above except `method` are required. The iterator can stream
frames without loading a complete trajectory into memory. Counts and
normalizations are pooled as described in
[Trajectory averaging](methods.md#trajectory-averaging).

## `plot_rdf`: plot a saved or computed curve

```python
from curdf import plot_rdf

fig = plot_rdf(
    bins,
    gr,
    path=None,
    show=False,
    xlabel="r (Å)",
    ylabel="g(r)",
    title=None,
)
```

`bins` and `gr` must have matching lengths. With a path, the function saves a
300 dpi figure and returns its Matplotlib `Figure`. It closes the figure after
saving. Create the destination directory yourself when plotting independently
of `rdf()`.

The plotting module selects the noninteractive `Agg` backend. Save figures with
`path` rather than relying on `show=True` to open a window.
