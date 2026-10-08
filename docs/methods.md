# Methods

## What the RDF measures

The radial distribution function `g(r)` compares the number of neighbours in a
spherical shell with the number expected from the selected groups' number
density. The radial coordinate uses the same length units as the input,
normally angstroms. `g(r)` is dimensionless.

cuRDF builds a neighbour list with NVIDIA ALCHEMI Toolkit-Ops, computes pair
distances with the returned periodic image shifts, and accumulates a radial
histogram with PyTorch. It then divides the counts by spherical shell volumes
and the population normalization below.

## Radial bins

The effective lower edge is:

```text
r0 = max(r_min, r_min_floor)
```

There are `nbins` equal-width bins between `r0` and `r_max`. A bin includes its
lower edge and excludes its upper edge. Distances equal to `r_max` and pairs
whose source and target atom indices are identical are excluded.

For bin `k`, the returned radial coordinate and shell volume are:

```text
r_k = (edge_k + edge_(k+1)) / 2
shell_volume_k = (4 pi / 3) * (edge_(k+1)^3 - edge_k^3)
```

The default `r_min=1.0` omits distances below 1 angstrom. Pass `r_min=0.0` to
include shorter distances. With the default floor, the first edge will be
`1e-6`, rather than exactly zero. Increasing `nbins` changes histogram
resolution, not the neighbour cutoff.

## Same-species normalization

For a group containing `N` atoms in a cell of volume `V`, the normalization is
`N * (N - 1) / V`. Each atom has `N - 1` possible partners because self-pairs
are excluded.

With `half_fill=True`, each unordered pair contributes once. The counts are
multiplied by two when forming the RDF:

```text
g_k = 2 * count_k / (shell_volume_k * N * (N - 1) / V)
```

With `half_fill=False`, both pair directions are counted and the factor of two
is omitted. The two conventions should give the same RDF with a complete
neighbour list.

This finite-population convention differs from using `N * N / V`. For the same
counts, the `N - 1` convention produces a value larger by `N / (N - 1)`. The
difference matters most for small groups. Report the normalization when
comparing with other packages.

## Cross-species normalization

For distinct groups A and B, with populations `N_A` and `N_B`, cuRDF counts
ordered A-to-B pairs and uses:

```text
g_AB,k = count_AB,k / (shell_volume_k * N_A * N_B / V)
```

Cross-species calculations disable `half_fill` automatically. There is no
extra factor of two. Use disjoint groups for cross-group calculations. The
normalization does not correct for overlapping, unequal selections.

For a same-species calculation through `rdf()`, supply `species_a` and omit
`species_b`. This also ensures that MDAnalysis uses one atom group for both
sides of the calculation.

## Trajectory averaging

For selected frames `f`, cuRDF sums pair counts and population normalizations
before dividing:

```text
Z_f = N_f * (N_f - 1) / V_f            # same-species
Z_f = N_A,f * N_B,f / V_f              # cross-species
g_k = pair_factor * sum_f(count_f,k) / (shell_volume_k * sum_f(Z_f))
```

Here `pair_factor` is two for unique same-species pairs and one otherwise.
Every selected frame contributes once. No weights based on frame time are
applied.

For fixed populations and volume, this equals an arithmetic mean of the
per-frame RDFs. When volume or populations change, it is an average weighted by
`Z_f`. Keep the species definitions and bin edges consistent when comparing
curves. cuRDF returns one pooled curve and does not calculate standard errors.

## Cells and periodic boundaries

Cell matrices have lattice vectors as rows. Orthorhombic and triclinic matrices
are forwarded to Toolkit-Ops. Periodic image shift vectors are converted to
Cartesian displacements using the cell matrix before distances are calculated.

The high-level `rdf()` default is fully periodic. It does not automatically
adopt ASE's stored periodicity flags. Pass `pbc=tuple(atoms.pbc)` when that is
the intended convention. The MDAnalysis adapter obtains each frame's cell from
its box dimensions.

The normalization always uses the full cell volume and three-dimensional
spherical shells. Setting an axis to nonperiodic changes the neighbour search,
but does not add a slab, surface, or finite-boundary correction.

For a conventional bulk RDF, use a cutoff below half the smallest
perpendicular cell width, so the spherical shells fit within the periodic
cell. For an orthorhombic box, these widths are its three side lengths. A
cutoff suitable for one frame may be too large for another if the cell shrinks.
The API does not enforce this cutoff choice or apply a shell-volume correction
for larger cutoffs.

## Numerical and performance controls

`method="cell_list"` reduces neighbour-search work for larger systems.
`method="naive"` considers pairs directly and can have less setup overhead for
small systems. Both feed the same RDF histogram and normalization. The cutoff,
density, GPU, and system size affect which is faster.

The default tensor dtype is float32. Distances near bin boundaries can move
between adjacent bins with changes in precision. `torch_dtype=torch.float64`
changes the tensor calculations, while the ASE and MDAnalysis adapters first
convert coordinates and cells to float32. Use `compute_rdf()` with float64
input arrays when preserving their input precision matters.

Increase `max_neighbors` if Toolkit-Ops reports insufficient neighbour capacity.
Use the same cutoff, bins, precision, and frame selection when comparing
methods or benchmarking. The package does not remove bonded pairs, exclude
same-molecule neighbours, or estimate coordination numbers automatically.
