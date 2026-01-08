import importlib
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pytest
import torch

# Ensure repo-local curdf is importable when running tests in-tree.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

rdf_module = importlib.import_module("curdf.rdf")


@pytest.fixture
def stub_neighbor(monkeypatch) -> List[Dict[str, object]]:
    """
    Replace Toolkit-Ops neighbor list with a simple CPU stub for tests.

    Records each call (method, half_fill, max_neighbors, cell, pbc) and returns
    deterministic neighbor pairs and zero shifts.
    """
    calls: List[Dict[str, object]] = []

    def _fake_neighbor_list(
        positions: torch.Tensor,
        r_max: float,
        cell: torch.Tensor,
        pbc: torch.Tensor,
        half_fill: bool = True,
        max_neighbors: int | None = None,
        method: str = "cell_list",
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        calls.append(
            {
                "method": method,
                "half_fill": half_fill,
                "max_neighbors": max_neighbors,
                "cell": cell.detach().cpu().numpy(),
                "pbc": pbc.detach().cpu().numpy(),
            }
        )
        pos = positions.detach().cpu().numpy()
        n_atoms = pos.shape[0]
        pairs: List[tuple[int, int]] = []
        if half_fill:
            for i in range(n_atoms):
                for j in range(i + 1, n_atoms):
                    if np.linalg.norm(pos[j] - pos[i]) <= r_max + 1e-6:
                        pairs.append((i, j))
        else:
            for i in range(n_atoms):
                for j in range(n_atoms):
                    if i == j:
                        continue
                    if np.linalg.norm(pos[j] - pos[i]) <= r_max + 1e-6:
                        pairs.append((i, j))

        if pairs:
            nlist = torch.tensor(pairs, device=positions.device, dtype=torch.int64).t()
        else:
            nlist = torch.empty((2, 0), device=positions.device, dtype=torch.int64)
        shifts = torch.zeros((nlist.shape[1], 3), device=positions.device, dtype=positions.dtype)
        return nlist, shifts

    monkeypatch.setattr(rdf_module, "build_neighbor_list", _fake_neighbor_list)
    return calls
