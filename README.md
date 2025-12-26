# cuRDF

GPU-accelerated radial distribution functions using Toolkit-Ops neighbor lists and PyTorch, with adapters for MDAnalysis and ASE.

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

ASE:
```python
from ase.io import read
from curdf import rdf_from_ase

atoms = read("POSCAR")
bins, gr = rdf_from_ase(atoms, selection=None, r_min=1.0, r_max=8.0, nbins=200)
```

## CLI
MDAnalysis:
```
rdf-gpu --format mdanalysis --topology top.data --trajectory traj.dcd --selection "name C" --r-max 8 --nbins 200 --device cuda --out results/rdf.npz --plot results/rdf.png
```

ASE:
```
rdf-gpu --format ase --ase-file POSCAR --selection 0,1,2 --r-max 8 --nbins 200 --device cuda
```

`--ordered-pairs` switches to counting ordered pairs (disable half-fill). `--no-wrap` leaves coordinates unwrapped if you already wrapped them upstream.

## Docs / examples / tests
- Docs in `docs/` (index, quickstart, api).
- Examples in `examples/` for basic, ASE, and MDAnalysis workflows.
- Tests in `tests/` (run with `pytest` or `pip install -e .[dev]` first).
