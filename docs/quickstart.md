# Quickstart

Install (editable for development):
```
pip install -e .[analysis]
```

Compute g(r) from an MDAnalysis trajectory:
```python
import MDAnalysis as mda
from curdf import rdf_from_mdanalysis

u = mda.Universe("top.data", "traj.dcd")
bins, gr = rdf_from_mdanalysis(u, selection="name C", r_min=1.0, r_max=8.0, nbins=200)
```

Compute g(r) from ASE:
```python
from ase.io import read
from curdf import rdf_from_ase

atoms = read("POSCAR")
bins, gr = rdf_from_ase(atoms, r_min=1.0, r_max=8.0, nbins=200)
```

CLI:
```
rdf-gpu --format mdanalysis --topology top.data --trajectory traj.dcd --selection "name C" --r-max 8 --nbins 200 --device cuda --plot rdf.png --out rdf.npz
```

For ASE input:
```
rdf-gpu --format ase --ase-file POSCAR --selection 0,1,2 --r-max 8 --nbins 200
```
