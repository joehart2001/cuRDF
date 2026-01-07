# API

## Python API

- `curdf.rdf(obj, species_a, species_b=None, neighbor_method="cell_list", **kwargs)`
  - Unified helper: pass an ASE Atoms/trajectory or an MDAnalysis Universe. `species_a` required; `species_b` defaults to same-species. Cross-species automatically disables half_fill. Common kwargs: `r_min`, `r_max`, `nbins`, `device`, `torch_dtype`, `max_neighbors`.
- `curdf.compute_rdf(...)`
  - Lower-level single-frame API (positions + cell + pbc), if you want to bypass adapters.
- `curdf.plot_rdf(bins, gr, path=None, show=False, xlabel="r (Å)", ylabel="g(r)", title=None)`
  - Quick plotting helper; saves if `path` provided, returns the Matplotlib figure.
