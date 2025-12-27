# Quickstart

Install (editable for development):
```
pip install -e .
```

Compute g(r) from ASE (XYZ/extxyz/ASE .traj):
```python
from ase.io import read
from curdf import rdf

atoms = read("structure.xyz")

# Unified helper: dispatches ASE vs MDAnalysis based on object type
bins, gr = rdf(atoms, species_a="C", species_b="O", r_min=1.0, r_max=8.0, nbins=200)
# Only use specific frames (e.g., first 1000 of a trajectory)
bins, gr = rdf(atoms_trajectory, species_a="C", species_b="O", r_min=1.0, r_max=8.0, nbins=200, index=":1000")
```

Compute g(r) from an MDAnalysis trajectory:
```python
import MDAnalysis as mda
from curdf import rdf

u = mda.Universe("top.data", "traj.dcd")
bins, gr = rdf(u, species_a="C", species_b="O", r_min=1.0, r_max=8.0, nbins=200)
```

CLI:
```
curdf --topology top.data --trajectory traj.dcd --species-a C --min 1 --max 8 --nbins 200 --device cuda --plot rdf.png --out rdf.npz
```

For ASE input (XYZ/extxyz/ASE .traj/LAMMPS data or dump via ASE readers):
```
curdf --file structure.xyz --species-a C --min 1 --max 8 --nbins 200
```

Cross-species via CLI: use `--species-a` and `--species-b` (or `--selection-a`/`--selection-b` for manual groups)
```
curdf --file structure.xyz --species-a C --species-b O --min 1 --max 8 --nbins 200 --device cuda
```
(`--species-b`/`--selection-b` automatically disables half-fill / pair doubling.)

LAMMPS dump (lammpstrj) without a separate topology, via MDAnalysis:
```
curdf --file dump.lammpstrj --species-a C --species-b O --min 1 --max 8 --nbins 200
```

If the LAMMPS data file needs a specific atom_style, pass `--atom-style "id type x y z"` (default is the same).
