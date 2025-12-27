import argparse
import sys
from pathlib import Path

import numpy as np


def _parse_args():
    p = argparse.ArgumentParser(description="cuRDF: GPU RDF using Toolkit-Ops + PyTorch")
    p.add_argument(
        "--format",
        choices=["mdanalysis", "ase", "lammps-dump"],
        required=True,
        help="Input backend",
    )
    p.add_argument("--topology", help="Topology file (MDAnalysis)")
    p.add_argument("--trajectory", nargs="+", help="Trajectory file(s) (MDAnalysis)")
    p.add_argument("--ase-file", help="Structure/trajectory file readable by ASE")
    p.add_argument("--ase-index", default=":", help="ASE index (default all frames)")
    p.add_argument("--selection", default=None, help="(Deprecated) alias for --selection-a")
    p.add_argument("--selection-a", default=None, help="MDAnalysis selection or ASE comma-separated indices for group A")
    p.add_argument("--selection-b", default=None, help="MDAnalysis selection or ASE comma-separated indices for group B")
    p.add_argument("--r-min", type=float, default=1.0)
    p.add_argument("--r-max", type=float, default=6.0)
    p.add_argument("--nbins", type=int, default=100)
    p.add_argument("--device", default="cuda")
    p.add_argument("--dtype", choices=["float32", "float64"], default="float32")
    p.add_argument("--half-fill", action="store_true", default=True, help="Use unique pairs (identical species)")
    p.add_argument("--ordered-pairs", action="store_true", help="Disable half-fill to count ordered pairs")
    p.add_argument("--max-neighbors", type=int, default=2048)
    p.add_argument("--no-wrap", action="store_true", help="Skip wrapping positions into the cell")
    p.add_argument("--plot", type=Path, help="Optional PNG plot output")
    p.add_argument("--out", type=Path, default=Path("rdf.npz"), help="NPZ output path")
    return p.parse_args()


def main():
    args = _parse_args()
    torch_dtype = {"float32": "float32", "float64": "float64"}[args.dtype]
    half_fill = False if args.ordered_pairs else args.half_fill
    # Cross-species (selection-b) should use ordered pairs
    if args.selection_b and half_fill:
        half_fill = False

    if args.format == "mdanalysis":
        if args.topology is None or args.trajectory is None:
            sys.exit("For mdanalysis format, provide --topology and --trajectory")
        try:
            import MDAnalysis as mda
        except ImportError:
            sys.exit("MDAnalysis is required for --format mdanalysis")

        u = mda.Universe(args.topology, *args.trajectory)
        selection_a = args.selection_a or args.selection
        selection_b = args.selection_b
        if selection_a is None:
            selection_a = "all"
        from .adapters import rdf_from_mdanalysis

        bins, gr = rdf_from_mdanalysis(
            u,
            selection=selection_a,
            selection_b=selection_b,
            r_min=args.r_min,
            r_max=args.r_max,
            nbins=args.nbins,
            device=args.device,
            torch_dtype=getattr(__import__("torch"), torch_dtype),
            half_fill=half_fill,
            max_neighbors=args.max_neighbors,
            wrap_positions=not args.no_wrap,
        )
    elif args.format == "lammps-dump":
        if args.trajectory is None:
            sys.exit("For lammps-dump format, provide --trajectory (LAMMPS dump / lammpstrj)")
        try:
            import MDAnalysis as mda
        except ImportError:
            sys.exit("MDAnalysis is required for --format lammps-dump")

        try:
            u = mda.Universe(args.trajectory[0], format="LAMMPSDUMP")
        except Exception as exc:
            sys.exit(f"Failed to load LAMMPS dump: {exc}")

        selection_a = args.selection_a or args.selection
        selection_b = args.selection_b
        if selection_a is None:
            selection_a = "all"
        from .adapters import rdf_from_mdanalysis

        bins, gr = rdf_from_mdanalysis(
            u,
            selection=selection_a,
            selection_b=selection_b,
            r_min=args.r_min,
            r_max=args.r_max,
            nbins=args.nbins,
            device=args.device,
            torch_dtype=getattr(__import__("torch"), torch_dtype),
            half_fill=half_fill,
            max_neighbors=args.max_neighbors,
            wrap_positions=not args.no_wrap,
        )
    else:
        if args.ase_file is None:
            sys.exit("For ase format, provide --ase-file")
        try:
            import ase.io
        except ImportError:
            sys.exit("ASE is required for --format ase")

        allowed_ext = {".xyz", ".extxyz", ".traj"}
        if Path(args.ase_file).suffix.lower() not in allowed_ext:
            sys.exit(f"ASE mode supports {sorted(allowed_ext)}; got {args.ase_file}")

        frames = ase.io.read(args.ase_file, index=args.ase_index)
        if isinstance(frames, list):
            atoms_or_traj = frames
        else:
            atoms_or_traj = frames

        sel_a = None
        sel_b = None
        selection_a = args.selection_a or args.selection
        if selection_a:
            sel_a = [int(x) for x in selection_a.split(",") if x.strip()]
        if args.selection_b:
            sel_b = [int(x) for x in args.selection_b.split(",") if x.strip()]

        from .adapters import rdf_from_ase

        bins, gr = rdf_from_ase(
            atoms_or_traj,
            selection=sel_a,
            selection_b=sel_b,
            r_min=args.r_min,
            r_max=args.r_max,
            nbins=args.nbins,
            device=args.device,
            torch_dtype=getattr(__import__("torch"), torch_dtype),
            half_fill=half_fill,
            max_neighbors=args.max_neighbors,
            wrap_positions=not args.no_wrap,
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.out, bins=bins, gr=gr)

    if args.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            sys.exit("matplotlib required for plotting")
        plt.plot(bins, gr)
        plt.xlabel("r (A)")
        plt.ylabel("g(r)")
        plt.hlines(1.0, xmin=args.r_min, xmax=args.r_max, colors="k", linestyles="dashed")
        plt.savefig(args.plot, dpi=300)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
