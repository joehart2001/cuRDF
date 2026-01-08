import math

import numpy as np
import torch

from curdf.rdf import accumulate_rdf, compute_rdf


def test_compute_rdf_single_frame_two_atoms(stub_neighbor):
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
    assert math.isclose(bins[1], 1.5)
    shell_vol = (4.0 / 3.0) * math.pi * (2.0**3 - 1.0**3)
    volume = 10.0**3
    n_atoms = 2
    rho = n_atoms / volume
    expected = (2.0 * 1.0) / (shell_vol * (n_atoms * rho))
    assert math.isclose(gr[1], expected, rel_tol=1e-5)


def test_accumulate_rdf_multiple_frames(stub_neighbor):
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
    assert gr[2] > 0
    assert gr[3] > 0


def test_compute_rdf_cross_species(stub_neighbor):
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
        half_fill=False,
        group_a_indices=[0],
        group_b_indices=[1],
    )
    shell_vol = (4.0 / 3.0) * math.pi * (2.0**3 - 1.0**3)
    volume = 10.0**3
    n_a = 1
    n_b = 1
    expected = 1.0 / (shell_vol * (n_a * (n_b / volume)))
    assert math.isclose(gr[1], expected, rel_tol=1e-5)


def test_compute_rdf_forwards_method(stub_neighbor):
    positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float32)
    cell = np.diag([6.0, 6.0, 6.0])
    compute_rdf(
        positions,
        cell,
        pbc=(True, True, True),
        r_min=0.0,
        r_max=3.0,
        nbins=3,
        device="cpu",
        torch_dtype=torch.float32,
        half_fill=True,
        method="naive",
    )
    assert stub_neighbor
    assert stub_neighbor[0]["method"] == "naive"


def test_accumulate_rdf_forwards_method_all_frames(stub_neighbor):
    cell = np.diag([5.0, 5.0, 5.0])
    frames = [
        {"positions": np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float32), "cell": cell, "pbc": (True, True, True)},
        {"positions": np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], dtype=np.float32), "cell": cell, "pbc": (True, True, True)},
    ]
    accumulate_rdf(
        frames,
        r_min=0.0,
        r_max=3.0,
        nbins=3,
        device="cpu",
        torch_dtype=torch.float32,
        half_fill=True,
        max_neighbors=128,
        method="naive",
    )
    assert len(stub_neighbor) == len(frames)
    assert all(call["method"] == "naive" for call in stub_neighbor)


def test_compute_rdf_triclinic_cell_forwarded(stub_neighbor):
    positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float32)
    cell = np.array(
        [
            [3.0, 0.1, 0.2],
            [0.0, 2.5, 0.3],
            [0.0, 0.0, 4.0],
        ],
        dtype=np.float32,
    )
    compute_rdf(
        positions,
        cell,
        pbc=(True, True, True),
        r_min=0.0,
        r_max=3.0,
        nbins=3,
        device="cpu",
        torch_dtype=torch.float32,
        half_fill=True,
        method="naive",
    )
    recorded_cell = stub_neighbor[0]["cell"]
    assert recorded_cell.shape[-2:] == (3, 3)
    assert np.allclose(recorded_cell[0], cell)


def test_compute_rdf_partial_pbc_forwarded(stub_neighbor):
    positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float32)
    cell = np.diag([6.0, 6.0, 6.0])
    compute_rdf(
        positions,
        cell,
        pbc=(True, False, True),
        r_min=0.0,
        r_max=3.0,
        nbins=3,
        device="cpu",
        torch_dtype=torch.float32,
        half_fill=True,
        method="naive",
    )
    recorded_pbc = stub_neighbor[0]["pbc"]
    assert recorded_pbc.shape[-1] == 3
    assert recorded_pbc[0, 0] and not recorded_pbc[0, 1] and recorded_pbc[0, 2]
