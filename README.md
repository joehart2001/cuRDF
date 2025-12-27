# cuRDF

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1085332119.svg)](https://doi.org/10.5281/zenodo.1085332119)


CUDA-accelerated radial distribution functions using NVIDIA ALCHEMI Toolkit-Ops O(N) neighbor lists and PyTorch. Compatible with ASE and MDAnalysis.

[![PyPI version](https://badge.fury.io/py/curdf.svg)](https://badge.fury.io/py/curdf)
[![Tests](https://github.com/josephhart/amorphous-carbon/actions/workflows/tests.yml/badge.svg)](https://github.com/josephhart/amorphous-carbon/actions/workflows/tests.yml)

## Install (editable)
```
pip install -e .
```
Add `[analysis]` extras if you want MDAnalysis/ASE/matplotlib:
```
pip install -e .[analysis]
```

## Library usage
```python
import curdf
import MDAnalysis as mda

u = mda.Universe("top.data", "traj.dcd")
bins, gr = curdf.rdf_from_mdanalysis(u, selection="name C", r_min=1.0, r_max=8.0, nbins=200)
```

ASE first (XYZ/extxyz/ASE .traj):
```python
from ase.io import read
from curdf import rdf_from_ase

atoms = read("structure.xyz")
bins, gr = rdf_from_ase(atoms, selection=None, r_min=1.0, r_max=8.0, nbins=200)  # selection=None -> all atoms
```

Cross-species (ASE): provide group A/B indices
```python
bins, gr = rdf_from_ase(atoms, selection=[0,1,2], selection_b=[3,4,5], r_min=1.0, r_max=8.0, nbins=200, half_fill=False)
```

MDAnalysis (explicit dependency required; also supports LAMMPS dump):
```python
import MDAnalysis as mda
from curdf import rdf_from_mdanalysis

u = mda.Universe("top.data", "traj.dcd")
bins, gr = curdf.rdf_from_mdanalysis(u, selection="name C", r_min=1.0, r_max=8.0, nbins=200)
```

## CLI
ASE (XYZ/extxyz/ASE .traj):
```
rdf-gpu --format ase --ase-file structure.xyz --selection 0,1,2 --r-max 8 --nbins 200 --device cuda
```

Cross-species via CLI (ASE indices or MDAnalysis selections):
```
rdf-gpu --format ase --ase-file structure.xyz --selection-a 0,1,2 --selection-b 3,4,5 --r-max 8 --nbins 200 --device cuda --ordered-pairs
```
(`--selection-b` automatically disables half-fill so pairs are ordered.)

LAMMPS dump (lammpstrj) via MDAnalysis:
```
rdf-gpu --format lammps-dump --trajectory dump.lammpstrj --selection "all" --r-max 8 --nbins 200 --device cuda
```

MDAnalysis:
```
rdf-gpu --format mdanalysis --topology top.data --trajectory traj.dcd --selection "name C" --r-max 8 --nbins 200 --device cuda --out results/rdf.npz --plot results/rdf.png
```

`--ordered-pairs` switches to counting ordered pairs (disable half-fill). `--no-wrap` leaves coordinates unwrapped if you already wrapped them upstream.

## Docs / examples / tests
- Docs in `docs/` (index, quickstart, api).
- Examples in `examples/` for basic, ASE, and MDAnalysis workflows.
- Tests in `tests/` (run with `pytest` or `pip install -e .[dev]` first).
- Build Sphinx docs with `pip install -e .[docs]` then `sphinx-build -b html docs/source docs/build/html` (footer: "Built with Sphinx using a theme provided by Read the Docs.").

## Citation
See `CITATION.cff` for how to cite cuRDF in your work.
