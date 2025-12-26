"""
Example: compute g(r) from an ASE Atoms object (in-memory).
Replace the generated Atoms with `ase.io.read("structure.xyz")` for real data (XYZ/extxyz).
"""

import numpy as np
from ase import Atoms
from ase.build import bulk

from curdf import rdf_from_ase


def main():
    # Simple silicon bulk (diamond) as a demo
    atoms = bulk("Si", "diamond", a=5.43).repeat((2, 2, 2))
    bins, gr = rdf_from_ase(
        atoms,
        selection=None,
        r_min=0.0,
        r_max=6.0,
        nbins=200,
        device="cpu",
        half_fill=True,
    )
    print("bins shape:", bins.shape, "g(r) shape:", gr.shape)
    print("first few g(r):", gr[:10])


if __name__ == "__main__":
    main()
