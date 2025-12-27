# Quickstart

Install (editable for development):
```
pip install -e .
```

Compute g(r) from ASE (XYZ/extxyz/ASE .traj):
```python
from ase.io import read
from curdf import rdf_from_ase

atoms = read("structure.xyz")
# selection=None means use all atoms; pass a list of indices to subset
bins, gr = rdf_from_ase(atoms, selection=None, r_min=1.0, r_max=8.0, nbins=200)
```

Cross-species RDF (ASE): provide group A/B indices
```python
bins, gr = rdf_from_ase(
    atoms,
    selection=[0, 1, 2],
    selection_b=[3, 4, 5],
    r_min=1.0,
    r_max=8.0,
    nbins=200,
    half_fill=False,  # cross-species -> ordered pairs
)
```

Compute g(r) from an MDAnalysis trajectory:
```python
import MDAnalysis as mda
from curdf import rdf_from_mdanalysis

u = mda.Universe("top.data", "traj.dcd")
bins, gr = rdf_from_mdanalysis(u, selection="name C", r_min=1.0, r_max=8.0, nbins=200)
```

CLI:
```
rdf-gpu --format mdanalysis --topology top.data --trajectory traj.dcd --selection "name C" --r-max 8 --nbins 200 --device cuda --plot rdf.png --out rdf.npz
```

For ASE input (XYZ/extxyz/ASE .traj):
```
rdf-gpu --format ase --ase-file structure.xyz --selection 0,1,2 --r-max 8 --nbins 200
```

Cross-species via CLI: use `--selection-a` and `--selection-b` (MDAnalysis selections or ASE index lists)
```
rdf-gpu --format ase --ase-file structure.xyz --selection-a 0,1,2 --selection-b 3,4,5 --r-max 8 --nbins 200 --device cuda --ordered-pairs
```
(`--selection-b` automatically disables half-fill / pair doubling.)

LAMMPS dump (lammpstrj) without a separate topology, via MDAnalysis:
```
rdf-gpu --format lammps-dump --trajectory dump.lammpstrj --selection "all" --r-max 8 --nbins 200
```
