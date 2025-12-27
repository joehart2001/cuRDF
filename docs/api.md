# API

## Python API

- `curdf.rdf(obj, species_a, species_b=None, **kwargs)`
  - Unified helper: pass an ASE Atoms/trajectory or an MDAnalysis Universe. `species_a` required; `species_b` defaults to same-species. Cross-species automatically disables half_fill. Common kwargs: `r_min`, `r_max`, `nbins`, `device`, `torch_dtype`, `max_neighbors`.
- `curdf.compute_rdf(...)`
  - Lower-level single-frame API (positions + cell + pbc), if you want to bypass adapters.

## CLI

`curdf --format mdanalysis --topology top.data --trajectory traj.dcd --species-a C --r-max 8 --nbins 200 --device cuda --out rdf.npz --plot rdf.png`

Key flags:
- `--format`: optional; auto-detected from inputs (topology+trajectory → mdanalysis; file extensions .xyz/.extxyz/.traj/.data/.dump → ase; .lammpstrj → lammps-dump)
- `--file`: ASE/LAMMPSDUMP file path (XYZ/extxyz/traj/LAMMPS data/dump)
- `--species-a` (required), `--species-b` (defaults to A): element symbols for groups A/B. `--selection-a`, `--selection-b`: manual selections/indices. `--selection` remains an alias for A.
- `--selection`: MDAnalysis selection string or comma-separated indices for ASE
- `--min/--max/--nbins`
- `--device`, `--dtype`, `--max-neighbors`
- `--no-wrap` if upstream wrapping already performed
