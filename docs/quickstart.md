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
# Same-species (specify both): ordered pairs collapsed internally
bins, gr = rdf_from_ase(atoms, species_a="C", species_b="C", r_min=1.0, r_max=8.0, nbins=200)

# Cross-species: specify element symbols; ordered pairs handled automatically
bins, gr = rdf_from_ase(atoms, species_a="C", species_b="O", r_min=1.0, r_max=8.0, nbins=200)
```

Compute g(r) from an MDAnalysis trajectory:
```python
import MDAnalysis as mda
from curdf import rdf_from_mdanalysis

u = mda.Universe("top.data", "traj.dcd")
bins, gr = rdf_from_mdanalysis(u, species_a="C", species_b="C", r_min=1.0, r_max=8.0, nbins=200)  # same-species
bins, gr = rdf_from_mdanalysis(u, species_a="C", species_b="O", r_min=1.0, r_max=8.0, nbins=200)  # cross-species
```

CLI:
```
curdf --format mdanalysis --topology top.data --trajectory traj.dcd --species-a C --r-max 8 --nbins 200 --device cuda --plot rdf.png --out rdf.npz
```

For ASE input (XYZ/extxyz/ASE .traj):
```
curdf --format ase --file structure.xyz --species-a C --r-max 8 --nbins 200
```

Cross-species via CLI: use `--species-a` and `--species-b` (or `--selection-a`/`--selection-b` for manual groups)
```
curdf --format ase --file structure.xyz --species-a C --species-b O --r-max 8 --nbins 200 --device cuda
```
(`--species-b`/`--selection-b` automatically disables half-fill / pair doubling.)

LAMMPS dump (lammpstrj) without a separate topology, via MDAnalysis:
```
curdf --format lammps-dump --file dump.lammpstrj --species-a C --species-b O --r-max 8 --nbins 200
```
