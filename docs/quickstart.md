# Quickstart

## Installation

Install the latest release from PyPI:

```bash
pip install cuRDF
```

For development, clone the repository and install it in editable mode:

```bash
git clone https://github.com/joehart2001/curdf.git
cd curdf
pip install -e .
```

Use Python 3.11 or newer and an NVIDIA GPU with a compatible CUDA installation.
The package is installed as `cuRDF` and imported as `curdf`. Check CUDA access
with:

```bash
python -c "import torch; print(torch.cuda.is_available())"
```

If this prints `False`, follow the [CUDA setup guidance](troubleshooting.md#cuda-is-not-available)
before running an RDF calculation.

## A single ASE structure

Read a structure with a cell and compute a cross-species RDF:

```python
from ase.io import read
from curdf import rdf, plot_rdf

atoms = read("structure.extxyz")
bins, gr = rdf(
    atoms,
    species_a="C",
    species_b="O",
    r_min=0.0,
    r_max=6.0,
    nbins=200,
    pbc=tuple(atoms.pbc),
    output="results/rdf_CO.npz",
)
plot_rdf(bins, gr, path="results/rdf_CO.png")
```

The structure must contain both selected species. `bins` contains radial bin
centres in angstroms and `gr` contains dimensionless values of `g(r)`. Both are
NumPy arrays of length `nbins` on the CPU, even though the calculation uses CUDA.

For a same-species RDF, omit `species_b`:

```python
bins, gr = rdf(atoms, species_a="C", r_min=0.0, r_max=6.0, nbins=200)
```

The default periodicity is `(True, True, True)`. Pass `pbc` explicitly to use
the flags stored by ASE. Choose a cutoff appropriate for the cell, as explained
in [Methods](methods.md#cells-and-periodic-boundaries).

## An ASE trajectory

Read the frames you want to analyse before passing them to `rdf()`:

```python
from ase.io import read
from curdf import rdf

# Skip the first 1000 frames and keep every tenth remaining frame.
frames = read("trajectory.extxyz", index="1000::10")
bins, gr = rdf(
    frames,
    species_a="C",
    r_min=0.0,
    r_max=6.0,
    nbins=200,
    output="results/rdf_CC.csv",
)
```

Pass a materialized list of frames. A filename by itself is not accepted by
`rdf()`. The returned curve pools counts across the selected frames, with the
normalization described in [Trajectory averaging](methods.md#trajectory-averaging).

## An MDAnalysis trajectory

Create a universe from the topology and trajectory, then select frames with
`start`, `stop`, and `step`:

```python
import MDAnalysis as mda
from curdf import rdf

u = mda.Universe("topology.data", "trajectory.dcd", atom_style="id type x y z")
bins, gr = rdf(
    u,
    species_a="O",
    species_b="H",
    atom_types_map={1: "O", 2: "H"},
    start=1000,
    step=10,
    r_min=0.0,
    r_max=6.0,
    nbins=200,
)
```

Adjust the type mapping to match your topology. MDAnalysis selections use atom
names. The mapping supplies names only when the topology has none. Existing
names such as `OW` and `HW` should instead be passed as `species_a="OW"` and
`species_b="HW"`. Each frame must have valid box dimensions.

## Saving and reusing results

`output` creates parent directories and writes the requested file. CSV and TSV
contain columns `r` and `g_r`. NPZ files contain arrays `bins` and `gr`:

```python
import numpy as np
from curdf import plot_rdf

data = np.load("results/rdf_CO.npz")
plot_rdf(data["bins"], data["gr"], path="results/rdf_CO_replot.png")
```

Use [API](api.md) for parameter defaults and other output formats, and
[Troubleshooting](troubleshooting.md) for setup or input errors.
