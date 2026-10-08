# Troubleshooting

## CUDA is not available

Check the Python environment that will run your calculation:

```bash
python -c "import torch; print(torch.__version__); print(torch.version.cuda); print(torch.cuda.is_available())"
nvidia-smi
```

cuRDF defaults to CUDA and requires an NVIDIA GPU for its production
neighbour-list backend. There is no automatic CPU fallback.
`torch.version.cuda` being `None` indicates a CPU-only PyTorch installation.
If CUDA is installed but unavailable, check the driver and GPU allocation
inside your job or container.

Install a PyTorch build appropriate for your machine using the
[official PyTorch selector](https://docs.pytorch.org/get-started/locally/).
Check the CUDA and dependency guidance in the
[Toolkit-Ops installation documentation](https://github.com/NVIDIA/nvalchemi-toolkit-ops).
Use Python 3.11 or newer for Toolkit-Ops. The exact driver and CUDA requirements
depend on the dependency versions you install.

## A species selection is empty

ASE and MDAnalysis use different labels. Inspect them before choosing species:

```python
# ASE
print(sorted(set(atoms.get_chemical_symbols())))

# MDAnalysis, if the topology has atom names
print(sorted(set(u.atoms.names)))
```

ASE normally selects symbols such as `C`, `O`, and `H`. MDAnalysis selects
names, which may instead be `OW`, `HW1`, or `HW2`. `atom_types_map` does not
replace existing MDAnalysis names.

If the MDAnalysis topology has types but no names, map every type used by the
system:

```python
bins, gr = rdf(u, species_a="O", species_b="H", atom_types_map={1: "O", 2: "H"})
```

For a topology without names or a mapping, a same-species call can label all
atoms as `species_a`. Supply a mapping for a mixture to avoid analysing all
atoms as one species. For same-species RDFs, omit `species_b`.

## A cell is missing or invalid

An XYZ file may not store a simulation cell. Use an extended XYZ, ASE trajectory,
or another format that preserves it, or set the correct cell before analysis.
The box volume enters the RDF normalization, so it must describe your system.

Inspect `atoms.cell.array` for ASE or `u.trajectory.ts.dimensions` for
MDAnalysis. MDAnalysis needs lengths and angles for each analysed frame.
`compute_rdf()` requires a `(3, 3)` cell and three Boolean periodicity flags.

For a slab or a nonperiodic structure, passing `pbc=(True, True, False)` changes
the neighbour search. It does not correct the bulk shell normalization for
interfaces or open boundaries. See [Methods](methods.md#cells-and-periodic-boundaries).

## No frames or no possible pairs

An empty frame range, empty group, or same-species group with fewer than two
atoms can give zero total normalization. Check that your frame range is
nonempty and that every selected group contains the intended atoms.

Use ASE reader strings in `ase.io.read()`, not in `rdf()`:

```python
from ase.io import read
from curdf import rdf

frames = read("trajectory.extxyz", index="1000::10")
bins, gr = rdf(frames, species_a="C")
```

For a list or an MDAnalysis trajectory, use a Python slice. To analyse a single
ASE frame from a list, pass `frames[k]` directly. Avoid one-shot generators at
the high-level entry point, and use `accumulate_rdf()` for streaming frame
dictionaries.

## Neighbour capacity or GPU memory is insufficient

If Toolkit-Ops reports insufficient neighbour capacity, increase
`max_neighbors` and rerun:

```python
bins, gr = rdf(atoms, species_a="C", max_neighbors=4096)
```

A larger cutoff or a denser system can require more neighbours and more memory.
Reduce the cutoff only if that still covers the distances needed for your
analysis. Compare `"cell_list"` and `"naive"` on a representative frame. Reducing
`nbins` changes histogram memory but does not reduce neighbour-list size.

## The RDF differs from another package

First align the selected atoms, frames, length units, periodicity, radial range,
and bin edges. Then compare the population normalization. cuRDF uses `N - 1`
possible partners for a same-species RDF, rather than `N`.

Also check whether the reference excludes bonded or same-molecule pairs.
cuRDF excludes self-pairs but does not apply those additional exclusions.
See [Methods](methods.md) for the formulas and trajectory weighting.

## A saved file or plot is missing

`rdf(output=...)` creates parent directories. Use a supported extension such as
`.npz` or `.csv`. If both `output` and `outdir` are set, `output` takes priority.
CSV and TSV export uses pandas. If that optional import fails, install pandas
in the active environment or choose NPZ or JSON.

When calling `plot_rdf()` independently, create the destination directory first.
The plotting module uses the noninteractive `Agg` backend, so save with `path`
instead of expecting `show=True` to open a window.

## What the CI tests cover

The hosted GitHub test jobs use a CPU neighbour-list stub to check RDF
normalization and adapter behavior. CUDA integration tests skip when a GPU is
unavailable. Passing hosted checks does not confirm GPU or driver compatibility
on your machine.

For an issue report, include the cuRDF, Toolkit-Ops, and PyTorch versions, the
CUDA diagnostic above, your cell and species labels, and a small input that
reproduces the problem.
