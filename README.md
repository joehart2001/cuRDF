# cuRDF

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1085332119.svg)](https://doi.org/10.5281/zenodo.1085332119) [![PyPI](https://img.shields.io/pypi/v/cuRDF.svg)](https://pypi.org/project/cuRDF/) [![Docs](https://readthedocs.org/projects/curdf/badge/?version=latest)](https://curdf.readthedocs.io/en/latest/)


CUDA-accelerated radial distribution functions using NVIDIA ALCHEMI Toolkit-Ops O(N) neighbor lists and PyTorch. Compatible with ASE (most common) and MDAnalysis.

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
```python
from ase.io import read
from curdf import rdf_from_ase

# ASE file formats: XYZ, extxyz, traj
atoms = read("structure.extxyz")

# Same-species: uses half_fill=True for the neighbour list to avoid double counting and save compute
bins, gr = rdf_from_ase(
  atoms, 
  species_a="C",
  species_b="C",
  r_min=1.0,
  r_max=8.0,
  nbins=200 # resolution of RDF
)

# Cross-species
bins, gr = rdf_from_ase(
  atoms, 
  species_a="C",
  species_b="O",
  r_min=1.0,
  r_max=8.0,
  nbins=200
)
```



MDAnalysis (also supports LAMMPS dump):
```python
import MDAnalysis as mda
from curdf import rdf_from_mdanalysis

u = mda.Universe("top.data", "traj.dcd")
bins, gr = curdf.rdf_from_mdanalysis(u, species_a="C", species_b="O", r_min=1.0, r_max=8.0, nbins=200)
```

## CLI
ASE (XYZ/extxyz/ASE .traj):
```
curdf --format ase --file structure.xyz --species-a C --r-max 8 --nbins 200 --device cuda
```

Cross-species via CLI:
```
curdf --format ase --file structure.xyz --species-a C --species-b O --r-max 8 --nbins 200 --device cuda
```

LAMMPS dump (lammpstrj) via MDAnalysis:
```
curdf --format lammps-dump --file dump.lammpstrj --species-a C --species-b O --r-max 8 --nbins 200 --device cuda
```

MDAnalysis:
```
curdf --format mdanalysis --topology top.data --trajectory traj.dcd --species-a C --r-max 8 --nbins 200 --device cuda --out results/rdf.npz --plot results/rdf.png
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
