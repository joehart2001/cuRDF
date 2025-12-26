# API

## Python API

- `curdf.compute_rdf(positions, cell, pbc, r_min=1.0, r_max=6.0, nbins=100, device="cuda", torch_dtype=torch.float32, half_fill=True, max_neighbors=2048)`
  - positions `(N,3)` array-like
  - cell `(3,3)` triclinic allowed; PBC flags `(3,)`
  - returns `(bins, gr)` as NumPy arrays for one frame

- `curdf.rdf_from_mdanalysis(universe, selection="all", ...)`
  - streams MDAnalysis trajectory, optional wrapping, selection string

- `curdf.rdf_from_ase(atoms_or_trajectory, selection=None, ...)`
  - ASE Atoms or iterable of Atoms; selection by indices

## CLI

`rdf-gpu --format mdanalysis --topology top.data --trajectory traj.dcd --selection "name C" --r-max 8 --nbins 200 --device cuda --out rdf.npz --plot rdf.png`

Key flags:
- `--format`: `mdanalysis` or `ase`
- `--selection`: MDAnalysis selection string or comma-separated indices for ASE
- `--r-min/--r-max/--nbins`
- `--device`, `--dtype`, `--max-neighbors`
- `--ordered-pairs` to disable half-fill normalization
- `--no-wrap` if upstream wrapping already performed
