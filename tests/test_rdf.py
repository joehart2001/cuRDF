import math
import numpy as np
import torch

from curdf.rdf import compute_rdf, accumulate_rdf


def test_compute_rdf_single_frame_two_atoms():
    positions = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], dtype=np.float32)
    cell = np.diag([10.0, 10.0, 10.0])
    bins, gr = compute_rdf(
        positions,
        cell,
        pbc=(True, True, True),
        r_min=0.0,
        r_max=5.0,
        nbins=5,
        device="cpu",
        torch_dtype=torch.float32,
        half_fill=True,
    )
    # Distance 1.5 should fall into bin index 1 (edges 0,1,2,...)
    assert math.isclose(bins[1], 1.5)
    # Expected g(r) from analytic normalization
    shell_vol = (4.0 / 3.0) * math.pi * (2.0**3 - 1.0**3)
    volume = 10.0**3
    n_atoms = 2
    rho = n_atoms / volume
    expected = (2.0 * 1.0) / (shell_vol * (n_atoms * rho))
    assert math.isclose(gr[1], expected, rel_tol=1e-5)


def test_accumulate_rdf_multiple_frames():
    cell = np.diag([8.0, 8.0, 8.0])
    positions1 = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], dtype=np.float32)
    positions2 = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]], dtype=np.float32)
    frames = [
        {"positions": positions1, "cell": cell, "pbc": (True, True, True)},
        {"positions": positions2, "cell": cell, "pbc": (True, True, True)},
    ]
    bins, gr = accumulate_rdf(
        frames,
        r_min=0.0,
        r_max=5.0,
        nbins=5,
        device="cpu",
        torch_dtype=torch.float32,
        half_fill=True,
        max_neighbors=1024,
    )
    # Distances fall into bins near 1 and 3; check both nonzero
    assert gr[1] > 0
    assert gr[3] > 0


def test_compute_rdf_cross_species():
    # A at origin, B at 1.5
    positions = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], dtype=np.float32)
    cell = np.diag([10.0, 10.0, 10.0])
    bins, gr = compute_rdf(
        positions,
        cell,
        pbc=(True, True, True),
        r_min=0.0,
        r_max=5.0,
        nbins=5,
        device="cpu",
        torch_dtype=torch.float32,
        half_fill=False,  # cross-species -> no pair doubling
        group_a_indices=[0],
        group_b_indices=[1],
    )
    shell_vol = (4.0 / 3.0) * math.pi * (2.0**3 - 1.0**3)
    volume = 10.0**3
    n_a = 1
    n_b = 1
    expected = 1.0 / (shell_vol * (n_a * (n_b / volume)))
    assert math.isclose(gr[1], expected, rel_tol=1e-5)
