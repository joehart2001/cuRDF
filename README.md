# cuRDF

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1085332119.svg)](https://doi.org/10.5281/zenodo.1085332119) [![PyPI](https://img.shields.io/pypi/v/cuRDF.svg)](https://pypi.org/project/cuRDF/)


CUDA-accelerated radial distribution functions using NVIDIA ALCHEMI Toolkit-Ops O(N) neighbor lists and PyTorch. Compatible with ASE (most common) and MDAnalysis.

Check out NVIDIA ALCHEMI Toolkit-Ops [GitHub](https://github.com/NVIDIA/nvalchemi-toolkit-ops)

## Install
Latest release:
```
pip install cuRDF
```
For development:
```
git clone https://github.com/joehart2001/curdf.git
cd curdf
pip install -e .
```

## Quickstart
Atoms object via ASE:
```python
from ase.io import read
from curdf import rdf

# ASE file formats: XYZ, extxyz, traj, LAMMPS data/dump (via ASE readers)
atoms = read("md_run.extxyz")

# Unified helper (dispatches ASE vs MDAnalysis based on object type)
bins, gr = rdf(atoms, species_a="C", species_b="O", r_min=1.0, r_max=8.0, nbins=200)

# Only use specific frames (e.g., first 1000 of a trajectory)
bins, gr = rdf(atoms_trajectory, species_a="C", species_b="O", r_min=1.0, r_max=8.0, nbins=200, index=":1000")
```



Topology and trajectory via MDAnalysis:
```python
import MDAnalysis as mda
from curdf import rdf

u = mda.Universe("topology.data", "traj.dcd")
bins, gr = rdf(u, species_a="C", species_b="O", r_min=1.0, r_max=8.0, nbins=200)
```

## CLI
ASE (XYZ/extxyz/ASE .traj/LAMMPS data or dump)
```
curdf --file structure.xyz --species-a C --species-b O --min 1 --max 8 --nbins 200 --device cuda
```

LAMMPS dump (lammpstrj) via MDAnalysis:
```
curdf --file dump.lammpstrj --species-a C --species-b O --min 1 --max 8 --nbins 200 --device cuda
```

MDAnalysis:
```
curdf --topology top.data --trajectory traj.dcd --species-a C --species-b C --min 1 --max 8 --nbins 200 --device cuda --out results/rdf.npz --plot results/rdf.png
# If the LAMMPS data file needs a specific atom_style, pass --atom-style "id type x y z" (default)
# If the LAMMPS data file lacks atom names, map types to elements: --atom-types "1:C,2:O"
```

`--no-wrap` leaves coordinates unwrapped if you already wrapped them upstream. Half-fill is chosen automatically based on species (same-species → half-fill).

## Citation
If you use cuRDF in your work, please cite:
```
@software{cuRDF,
  author    = {Hart, Joseph},
  title     = {cuRDF: GPU-accelerated radial distribution functions},
  month     = dec,
  year      = 2025,
  publisher = {Zenodo},
  version   = {0.1.0},
  doi       = {10.5281/zenodo.1085332119},
  url       = {https://doi.org/10.5281/zenodo.1085332119}
}
```
