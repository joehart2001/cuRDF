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

# CLI has been removed; use the Python API shown above.
If the LAMMPS data file needs a specific atom_style, pass `--atom-style "id type x y z"` (default is the same).
