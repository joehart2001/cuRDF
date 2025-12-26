"""
Minimal in-memory RDF on CPU for quick verification.
"""

import numpy as np
import torch
from curdf import compute_rdf


def main():
    positions = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], dtype=np.float32)
    cell = np.diag([10.0, 10.0, 10.0]).astype(np.float32)
    bins, gr = compute_rdf(
        positions,
        cell,
        pbc=(True, True, True),
        r_min=0.0,
        r_max=5.0,
        nbins=50,
        device="cpu",
        torch_dtype=torch.float32,
        half_fill=True,
    )
    print("bins:", bins)
    print("g(r):", gr)


if __name__ == "__main__":
    main()
