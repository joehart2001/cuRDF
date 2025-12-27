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

atoms = read("structure.xyz")
bins, gr = rdf_from_ase(atoms, selection=None, r_min=1.0, r_max=8.0, nbins=200)  # selection=None -> all atoms
```

Cross-species (ASE): provide group A/B indices
```python
bins, gr = rdf_from_ase(
    atoms,
    selection=[0,1,2],   # group A
    selection_b=[3,4,5], # group B
    r_min=1.0,
    r_max=8.0,
    nbins=200,
    half_fill=False,     # ordered pairs for cross-species
)
```

ASE file formats: XYZ, extxyz, ASE .traj. (Other formats supported by ASE may work, but these are tested.)

MDAnalysis (also supports LAMMPS dump):
```python
import MDAnalysis as mda
from curdf import rdf_from_mdanalysis

u = mda.Universe("top.data", "traj.dcd")
bins, gr = curdf.rdf_from_mdanalysis(u, selection="name C", r_min=1.0, r_max=8.0, nbins=200)
```

Same vs cross-species (why it matters):
- Same-species default (`selection=None` or single selection) uses unique pairs only (`half_fill=True`), cutting neighbor storage and skipping duplicate A–B/B–A work.
- Cross-species: pass both groups (`selection` and `selection_b`); we use ordered pairs (`half_fill=False`) so only A–B contributes.
- Specifying species controls which pairs are built and how normalization is done. Same-species mode saves compute; cross-species needs explicit groups to target only those pairs.

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
- Hosted docs (Read the Docs): https://curdf.readthedocs.io/en/latest/
- Sources in `docs/` (index, quickstart, api).
- Examples in `examples/` for basic, ASE, and MDAnalysis workflows.
- Tests in `tests/` (run with `pytest`).
- Docs deploy via GitHub Pages (workflow `.github/workflows/docs.yml`) and are compatible with Read the Docs (`.readthedocs.yaml`); local build: `sphinx-build -b html docs/source docs/build/html`.

## Citation
DOI: https://doi.org/10.5281/zenodo.1085332119  
See `CITATION.cff` for how to cite cuRDF in your work.
