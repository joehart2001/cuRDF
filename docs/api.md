# API

## Python API (common cases first)

- `curdf.rdf_from_ase(atoms_or_trajectory, selection=None, selection_b=None, species_a, species_b=None, ...)`
  - ASE Atoms or iterable of Atoms. `species_a` is required; `species_b` defaults to same-species. `selection`/`selection_b` allow manual index lists. Cross-species automatically disables half_fill. Intended for XYZ/extxyz/ASE .traj inputs.

- `curdf.rdf_from_mdanalysis(universe, selection=None, selection_b=None, species_a, species_b=None, ...)`
  - Streams MDAnalysis trajectory, optional wrapping. `species_a` is required; `species_b` defaults to same-species. Cross-species automatically disables half_fill.

- `curdf.compute_rdf(positions, cell, pbc, r_min=1.0, r_max=6.0, nbins=100, device="cuda", torch_dtype=torch.float32, half_fill=True, max_neighbors=2048)`
  - Lower-level API: positions `(N,3)` array-like, cell `(3,3)` triclinic allowed; PBC flags `(3,)`; returns `(bins, gr)` for one frame.

## CLI

`curdf --format mdanalysis --topology top.data --trajectory traj.dcd --species-a C --r-max 8 --nbins 200 --device cuda --out rdf.npz --plot rdf.png`

Key flags:
- `--format`: `mdanalysis`, `lammps-dump` (MDAnalysis LAMMPSDUMP), or `ase` (XYZ/extxyz/.traj)
- `--file`: ASE/LAMMPSDUMP file path
- `--species-a` (required), `--species-b` (defaults to A): element symbols for groups A/B. `--selection-a`, `--selection-b`: manual selections/indices. `--selection` remains an alias for A.
- `--selection`: MDAnalysis selection string or comma-separated indices for ASE
- `--r-min/--r-max/--nbins`
- `--device`, `--dtype`, `--max-neighbors`
- `--no-wrap` if upstream wrapping already performed
