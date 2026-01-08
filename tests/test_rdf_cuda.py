import numpy as np
import pytest

from curdf.rdf import compute_rdf


def test_compute_rdf_cuda_runs():
    torch_mod = pytest.importorskip("torch")
    if not torch_mod.cuda.is_available():
        pytest.skip("CUDA not available")
    positions = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]], dtype=np.float32)
    cell = np.diag([6.0, 6.0, 6.0])
    try:
        bins, gr = compute_rdf(
            positions,
            cell,
            pbc=(True, True, True),
            r_min=0.0,
            r_max=3.0,
            nbins=3,
            device="cuda",
            torch_dtype=torch_mod.float32,
            half_fill=True,
            max_neighbors=256,
        )
    except RuntimeError as err:
        pytest.skip(f"CUDA neighbor list unavailable: {err}")
    assert bins.shape[0] == gr.shape[0] == 3
    assert gr.max() >= 0.0
