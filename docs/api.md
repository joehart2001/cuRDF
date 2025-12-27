# API

## Python API (common cases first)

- `curdf.rdf_from_ase(atoms_or_trajectory, selection=None, selection_b=None, ...)`
  - ASE Atoms or iterable of Atoms; selection by indices (None = all atoms). Provide `selection_b` for cross-species RDF (half_fill automatically disabled). Intended for XYZ/extxyz/ASE .traj inputs.

- `curdf.rdf_from_mdanalysis(universe, selection="all", selection_b=None, ...)`
  - Streams MDAnalysis trajectory, optional wrapping. Provide `selection_b` for cross-species RDF (half_fill automatically disabled).

- `curdf.compute_rdf(positions, cell, pbc, r_min=1.0, r_max=6.0, nbins=100, device="cuda", torch_dtype=torch.float32, half_fill=True, max_neighbors=2048)`
  - Lower-level API: positions `(N,3)` array-like, cell `(3,3)` triclinic allowed; PBC flags `(3,)`; returns `(bins, gr)` for one frame.

## CLI

`rdf-gpu --format mdanalysis --topology top.data --trajectory traj.dcd --selection "name C" --r-max 8 --nbins 200 --device cuda --out rdf.npz --plot rdf.png`

Key flags:
- `--format`: `mdanalysis`, `lammps-dump` (MDAnalysis LAMMPSDUMP), or `ase` (XYZ/extxyz/.traj)
- `--selection-a`, `--selection-b`: group A/B (MDAnalysis selections or ASE index lists). `--selection` remains an alias for A.
- `--selection`: MDAnalysis selection string or comma-separated indices for ASE
- `--r-min/--r-max/--nbins`
- `--device`, `--dtype`, `--max-neighbors`
- `--ordered-pairs` to disable half-fill normalization
- `--no-wrap` if upstream wrapping already performed
