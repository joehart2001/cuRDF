import numpy as np
import pytest
import torch

from curdf.rdf import rdf as unified_rdf


def test_rdf_with_ase_input_cross_species_forwards_method_and_half_fill(stub_neighbor):
    ase = pytest.importorskip("ase")
    from ase import Atoms

    atoms = Atoms(
        symbols=["H", "He"],
        positions=[[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]],
        cell=[8.0, 8.0, 8.0],
        pbc=True,
    )
    unified_rdf(
        atoms,
        species_a="H",
        species_b="He",
        r_min=0.0,
        r_max=3.0,
        nbins=4,
        device="cpu",
        torch_dtype=torch.float32,
        half_fill=True,
        max_neighbors=64,
        method="naive",
    )
    assert stub_neighbor
    assert stub_neighbor[0]["method"] == "naive"
    assert stub_neighbor[0]["half_fill"] is False


def test_rdf_with_mdanalysis_input_cross_species_forwards_method(stub_neighbor):
    mda = pytest.importorskip("MDAnalysis")
    from MDAnalysis.coordinates.memory import MemoryReader

    coords = np.array([[[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]], dtype=np.float32)
    reader = MemoryReader(
        coords,
        order="fac",
        dimensions=np.array([[8.0, 8.0, 8.0, 90.0, 90.0, 90.0]], dtype=np.float32),
    )
    universe = mda.Universe.empty(2)
    universe.add_TopologyAttr("name", ["A", "B"])
    universe.trajectory = reader

    unified_rdf(
        universe,
        species_a="A",
        species_b="B",
        r_min=0.0,
        r_max=3.0,
        nbins=4,
        device="cpu",
        torch_dtype=torch.float32,
        half_fill=True,
        max_neighbors=64,
        method="naive",
        wrap_positions=False,
    )
    assert stub_neighbor
    assert stub_neighbor[0]["method"] == "naive"
    assert stub_neighbor[0]["half_fill"] is False


def test_rdf_accepts_list_of_ase_atoms(stub_neighbor):
    ase = pytest.importorskip("ase")
    from ase import Atoms

    frame = Atoms(
        symbols=["H", "H"],
        positions=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        cell=[6.0, 6.0, 6.0],
        pbc=True,
    )
    traj = [frame, frame.copy()]
    unified_rdf(
        traj,
        species_a="H",
        r_min=0.0,
        r_max=3.0,
        nbins=3,
        device="cpu",
        torch_dtype=torch.float32,
        half_fill=True,
        max_neighbors=64,
        method="naive",
    )
    assert len(stub_neighbor) == len(traj)
    assert all(call["method"] == "naive" for call in stub_neighbor)
