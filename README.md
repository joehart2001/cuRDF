# cuRDF

[![DOI](https://zenodo.org/badge/1085332119.svg)](https://doi.org/10.5281/zenodo.1085332119) [![PyPI](https://img.shields.io/pypi/v/cuRDF.svg)](https://pypi.org/project/cuRDF/)


CUDA-accelerated radial distribution functions using [NVIDIA ALCHEMI Toolkit-Ops](https://github.com/NVIDIA/nvalchemi-toolkit-ops) O(N) neighbor lists and PyTorch. Compatible with ASE Atoms or MDAnalysis Universe objects.

## Benchmarking

![cuRDF benchmark](benchmarks/results/results.png)

cuRDF is benchmarked against RDF (MDAnalysis) and neighbour list implementations on CPU (AMD Ryzen 9 9950X, 32 threads) and GPU (NVIDIA RTX 5090) for systems of varying sizes at a density of 0.05 atoms/Å³ over 1000 frames.

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
ASE Atoms object:
```python
from ase.io import read
from curdf import rdf

# Load trajectory or frame e.g. XYZ, extxyz, traj, LAMMPS data/dump
atoms = read("md_run.extxyz")

# Compute RDF between species C and O from 1.0 to 8.0 Å
bins, gr = rdf(
  atoms,
  species_a="C",
  species_b="O", # species b can be the same as species a
  r_min=1.0,
  r_max=8.0,
  nbins=200, # resolution of rdf histogram binning
  method="cell_list", # neighbor list method: "cell_list" (larger systems) or "naive" (smaller systems, less overhead)
  output = "results/rdf.csv" # optional output
)

# Plot RDF
plot_rdf(bins, gr, path="results/rdf.png")
```



MDAnalysis Universe (topology and trajectory):
```python
import MDAnalysis as mda
from curdf import rdf

u = mda.Universe("topology.data", "traj.dcd", atom_style="id type x y z")
bins, gr = rdf(
  u,
  species_a="C",
  species_b="O",
  r_min=1.0,
  r_max=8.0,
)
```

If the topology lacks atom names (only numeric types), supply a mapping:
```python
bins, gr = rdf(u, species_a="C", species_b="O", atom_types_map={1: "C", 2: "H"})
```

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
