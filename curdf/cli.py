import argparse
import sys
from pathlib import Path

import numpy as np


def _parse_args():
    p = argparse.ArgumentParser(description="cuRDF: GPU RDF using Toolkit-Ops + PyTorch")
    p.add_argument(
        "--format",
        choices=["mdanalysis", "ase", "lammps-dump"],
        required=False,
        help="Input backend (optional; auto-detected from inputs)",
    )
    p.add_argument("--topology", help="Topology file (MDAnalysis)")
    p.add_argument("--trajectory", nargs="+", help="Trajectory file(s) (MDAnalysis)")
    p.add_argument("--file", help="Structure/trajectory file (ASE or LAMMPSDUMP)")
    p.add_argument("--ase-file", help="(Deprecated) use --file for ASE")
    p.add_argument("--ase-index", default=":", help="ASE index (default all frames)")
    p.add_argument("--selection", default=None, help="(Deprecated) alias for --selection-a")
    p.add_argument("--selection-a", default=None, help="MDAnalysis selection or ASE comma-separated indices for group A")
    p.add_argument("--selection-b", default=None, help="MDAnalysis selection or ASE comma-separated indices for group B")
    p.add_argument("--species-a", required=True, help="Element symbol for group A (required)")
    p.add_argument("--species-b", default=None, help="Element symbol for group B (defaults to group A)")
    p.add_argument("--min", dest="r_min", type=float, default=1.0, help="Minimum r")
    p.add_argument("--max", dest="r_max", type=float, default=6.0, help="Maximum r")
    p.add_argument("--r-min", dest="r_min", type=float, help="(Deprecated) use --min")
    p.add_argument("--r-max", dest="r_max", type=float, help="(Deprecated) use --max")
    p.add_argument("--nbins", type=int, default=100)
    p.add_argument("--device", default="cuda", help="Device string (e.g., cuda or cpu)")
    p.add_argument("--dtype", choices=["float32", "float64"], default="float32")
    p.add_argument("--half-fill", action="store_true", default=True, help="Use unique pairs (identical species); auto-set based on species")
    p.add_argument("--max-neighbors", type=int, default=2048)
    p.add_argument("--no-wrap", action="store_true", help="Skip wrapping positions into the cell")
    p.add_argument("--plot", type=Path, help="Optional PNG plot output")
    p.add_argument("--out", type=Path, default=Path("rdf.npz"), help="NPZ output path")
    return p.parse_args()


def main():
    args = _parse_args()
    torch_dtype = {"float32": "float32", "float64": "float64"}[args.dtype]
    half_fill = args.half_fill
    # Cross-species should use ordered pairs
    if (args.species_b and args.species_b != args.species_a) or args.selection_b:
        half_fill = False

    # Infer format if not provided
    fmt = args.format
    if fmt is None:
        if args.topology and args.trajectory:
            fmt = "mdanalysis"
        elif args.file:
            ext = Path(args.file).suffix.lower()
            if ext in {".xyz", ".extxyz", ".traj"}:
                fmt = "ase"
            elif ext in {".lammpstrj", ".dump"}:
                fmt = "lammps-dump"
            else:
                sys.exit(f"Could not infer format from file extension '{ext}'. Provide --format.")
        else:
            sys.exit("Provide --format or inputs to infer it (topology+trajectory, or --file).")

    if fmt == "mdanalysis":
        if args.topology is None or args.trajectory is None:
            sys.exit("For mdanalysis format, provide --topology and --trajectory")
        try:
            import MDAnalysis as mda
        except ImportError:
            sys.exit("MDAnalysis is required for --format mdanalysis")

        u = mda.Universe(args.topology, *args.trajectory)
        selection_a = args.selection_a or args.selection
        selection_b = args.selection_b
        from .adapters import rdf_from_mdanalysis

        bins, gr = rdf_from_mdanalysis(
            u,
            species_a=args.species_a,
            species_b=args.species_b,
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
    elif fmt == "lammps-dump":
        if args.trajectory is None:
            if args.file is None:
                sys.exit("For lammps-dump format, provide --file or --trajectory (LAMMPS dump / lammpstrj)")
            args.trajectory = [args.file]
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
        from .adapters import rdf_from_mdanalysis

        bins, gr = rdf_from_mdanalysis(
            u,
            species_a=args.species_a,
            species_b=args.species_b,
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
    else:  # ase
        target_file = args.file or args.ase_file
        if target_file is None:
            sys.exit("For ase format, provide --file")
        try:
            import ase.io
        except ImportError:
            sys.exit("ASE is required for --format ase")

        allowed_ext = {".xyz", ".extxyz", ".traj"}
        if Path(target_file).suffix.lower() not in allowed_ext:
            sys.exit(f"ASE mode supports {sorted(allowed_ext)}; got {target_file}")

        frames = ase.io.read(target_file, index=args.ase_index)
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
            species_a=args.species_a,
            species_b=args.species_b,
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
