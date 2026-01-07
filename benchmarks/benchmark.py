import argparse
import csv
import time
from pathlib import Path

import numpy as np
import torch
from ase import Atoms
from ase.io import write

from curdf import compute_rdf

try:
    from plotting import plot_elapsed
    from alt_methods import bench_mdanalysis, bench_matscipy, bench_ase, bench_pymatgen
except ImportError:
    from .plotting import plot_elapsed
    from .alt_methods import bench_mdanalysis, bench_matscipy, bench_ase, bench_pymatgen

SIZES = [100, 250, 500, 750, 1000, 2000, 3000, 4000, 5000, 7500, 10000, 20000]  # atom counts
#SIZES = [100, 500, 1000]
N_FRAMES = 1000
DENSITY = 0.05  # atoms per Å^3 for synthetic boxes
R_MIN = 1.0
R_MAX = 6.0
NBINS = 100


def make_frames(n_atoms: int, n_frames: int, density: float, seed: int = 0):
    """Generate random boxes at a fixed number density so scaling is controlled."""
    box_len = float((n_atoms / density) ** (1.0 / 3.0))
    rng = np.random.default_rng(seed)
    frames = [rng.random((n_atoms, 3), dtype=np.float32) * box_len for _ in range(n_frames)]
    cell = np.diag([box_len, box_len, box_len]).astype(np.float32)
    return frames, cell


def _sync_if_needed(device: str):
    """Force CUDA to finish work so wall-clock timing is accurate."""
    if device == "cuda":
        torch.cuda.synchronize()


def bench_device(device: str, dtype: torch.dtype, output_rows: list, warmup_iters: int):
    if device == "cuda" and not torch.cuda.is_available():
        print("Skipping CUDA: torch.cuda.is_available() is False")
        return

    print(f"\nBenchmarking device={device}")
    edges = np.linspace(R_MIN, R_MAX, NBINS + 1)
    for n_atoms in SIZES:
        frames, cell = make_frames(n_atoms, N_FRAMES, DENSITY, seed=42)

        methods: list[tuple[str, str]] = [("cuRDF (cell list)", "cell_list"), ("cuRDF (naive)", "naive")]
        if device == "cpu":
            methods.extend(
                [
                    ("MDAnalysis", "mdanalysis"),
                    ("matscipy", "matscipy"),
                    ("ASE", "ase"),
                    ("pymatgen", "pymatgen"),
                ]
            )

        for label, method_key in methods:
            if method_key == "mdanalysis":
                for _ in range(warmup_iters):
                    bench_mdanalysis(frames, cell, R_MIN, R_MAX, NBINS)
                elapsed = bench_mdanalysis(frames, cell, R_MIN, R_MAX, NBINS)
                if elapsed is None:
                    continue
            elif method_key == "matscipy":
                for _ in range(warmup_iters):
                    bench_matscipy(frames[:1], cell, edges, R_MIN, R_MAX)
                elapsed = bench_matscipy(frames, cell, edges, R_MIN, R_MAX)
                if elapsed is None:
                    continue
            elif method_key == "ase":
                for _ in range(warmup_iters):
                    bench_ase(frames[:1], cell, edges, R_MIN, R_MAX)
                elapsed = bench_ase(frames, cell, edges, R_MIN, R_MAX)
                if elapsed is None:
                    continue
            elif method_key == "pymatgen":
                for _ in range(warmup_iters):
                    bench_pymatgen(frames[:1], cell, edges, R_MIN, R_MAX)
                elapsed = bench_pymatgen(frames, cell, edges, R_MIN, R_MAX)
                if elapsed is None:
                    continue
            else:
                for _ in range(warmup_iters):
                    _ = compute_rdf(
                        frames[0],
                        cell,
                        pbc=(True, True, True),
                        r_min=R_MIN,
                        r_max=R_MAX,
                        nbins=NBINS,
                        device=device,
                        torch_dtype=dtype,
                        method=method_key,
                    )
                _sync_if_needed(device)
                start = time.perf_counter()
                for pos in frames:
                    _ = compute_rdf(
                        pos,
                        cell,
                        pbc=(True, True, True),
                        r_min=R_MIN,
                        r_max=R_MAX,
                        nbins=NBINS,
                        device=device,
                        torch_dtype=dtype,
                        method=method_key,
                    )
                _sync_if_needed(device)
                elapsed = time.perf_counter() - start

            fps = N_FRAMES / elapsed
            print(f"N={n_atoms:6d} | method={label:14s} | frames/s={fps:8.2f} | elapsed={elapsed:6.2f}s")
            output_rows.append(
                {
                    "device": device,
                    "dtype": str(dtype).replace("torch.", ""),
                    "method": label,
                    "n_atoms": n_atoms,
                    "n_frames": N_FRAMES,
                    "elapsed_s": elapsed,
                    "frames_per_s": fps,
                }
            )


def main():
    """Run benchmarks, write CSV + plots, and dump example XYZ frames."""
    parser = argparse.ArgumentParser(description="cuRDF synthetic benchmark")
    parser.add_argument("--out", type=Path, default=Path("results/results.csv"), help="CSV output path")
    parser.add_argument("--dtype", choices=["float32", "float64"], default="float32")
    parser.add_argument("--warmup", type=int, default=1, help="Warm-up iterations per (N, method)")
    args = parser.parse_args()

    dtype = getattr(torch, args.dtype)
    rows = []
    overall_start = time.perf_counter()
    bench_device("cpu", dtype=dtype, output_rows=rows, warmup_iters=args.warmup)
    bench_device("cuda", dtype=dtype, output_rows=rows, warmup_iters=args.warmup)
    total_elapsed = time.perf_counter() - overall_start

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["device", "dtype", "method", "n_atoms", "n_frames", "elapsed_s", "frames_per_s"],
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nSaved results to {args.out}")
    print(f"Total benchmark wall time: {total_elapsed / 60.0:0.2f} min")

    plot_path = args.out.with_suffix(".png")
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    plot_elapsed(args.out, plot_path, n_frames=N_FRAMES)

    # Save first frame of each system as XYZ in a subfolder
    frames_dir = args.out.parent / "structure_example_frames"
    frames_dir.mkdir(parents=True, exist_ok=True)
    for n_atoms in SIZES:
        frames, cell = make_frames(n_atoms, 1, DENSITY, seed=42)
        pos = frames[0]
        symbols = ["X"] * n_atoms
        atoms = Atoms(symbols=symbols, positions=pos, cell=cell, pbc=True)
        xyz_out = frames_dir / f"frame_{n_atoms}.xyz"
        write(xyz_out, atoms)
    print(f"Saved example frames to {frames_dir}")


if __name__ == "__main__":
    main()
