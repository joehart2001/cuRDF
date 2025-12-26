"""
Example: compute g(r) from MDAnalysis trajectory.
Replace topology/trajectory with your files.
"""

import sys

try:
    import MDAnalysis as mda
except ImportError:
    sys.exit("Install MDAnalysis to run this example: pip install mdanalysis")

from curdf import rdf_from_mdanalysis


def main():
    topo = "topology.data"  # replace
    traj = "trajectory.dcd"  # replace

    u = mda.Universe(topo, traj)
    bins, gr = rdf_from_mdanalysis(
        u,
        selection="name C",
        r_min=1.0,
        r_max=8.0,
        nbins=200,
        device="cuda",
        half_fill=True,
    )
    print("bins shape:", bins.shape, "g(r) shape:", gr.shape)


if __name__ == "__main__":
    main()
