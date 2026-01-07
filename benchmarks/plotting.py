import csv
from pathlib import Path
import argparse

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def plot_elapsed(csv_path: Path, plot_path: Path, n_frames: int):
    if not csv_path.exists():
        print(f"Cannot plot: CSV not found at {csv_path}")
        return

    series: dict[tuple[str, str], list[tuple[int, float]]] = {}
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        missing_cols = {"device", "n_atoms", "elapsed_s"} - set(reader.fieldnames or [])
        if missing_cols:
            print(f"Cannot plot: CSV missing columns {sorted(missing_cols)}")
            return
        label_map = {
            "cell_list": "cuRDF (cell list)",
            "naive": "cuRDF (naive)",
            "mdanalysis": "MDAnalysis",
            "matscipy": "MatSciPy",
            "ase": "ASE",
            "pymatgen": "Pymatgen",
        }
        for row in reader:
            device = row["device"]
            method_raw = row.get("method", "unspecified")
            method = label_map.get(method_raw, method_raw)
            n_atoms = int(float(row["n_atoms"]))
            elapsed = float(row["elapsed_s"])
            key = (device, method)
            series.setdefault(key, []).append((n_atoms, elapsed))

    if not series:
        print("No data rows found; skipping plot")
        return

    # Color keyed by method so CPU/GPU share color; linestyle encodes device.
    color_map: dict[str, str] = {}
    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    def _get_color(method: str):
        if method not in color_map:
            color_map[method] = color_cycle[len(color_map) % len(color_cycle)]
        return color_map[method]

    fig, (ax_full, ax_zoom) = plt.subplots(1, 2, figsize=(12, 5))

    def _sort_key(item):
        (device, method), _ = item
        priority = {"cuRDF (cell list)": 0, "cuRDF (naive)": 1}
        return (priority.get(method, 2), method, device)

    series_handles = {}
    for key, pts in sorted(series.items(), key=_sort_key):
        device, method = key
        pts_sorted = sorted(pts, key=lambda x: x[0])
        n_atoms, elapsed = zip(*pts_sorted)
        color = _get_color(method)
        linestyle = "--" if device == "cpu" else "-"
        marker = "s" if device == "cpu" else "o"
        label = method

        line_full = ax_full.plot(n_atoms, elapsed, marker=marker, linestyle=linestyle, color=color, label=label)[0]
        series_handles[label] = line_full

        n_zoom = [n for n in n_atoms if n <= 1000]
        e_zoom = [e for n, e in zip(n_atoms, elapsed) if n <= 1000]
        if n_zoom:
            ax_zoom.plot(n_zoom, e_zoom, marker=marker, linestyle=linestyle, color=color, label=label)

    for ax, title in ((ax_full, "All sizes"), (ax_zoom, "n_atoms < 1000")):
        ax.set_xlabel("Number of atoms")
        ax.set_ylabel(f"Elapsed time (s) for {n_frames} frames (log scale)")
        ax.set_yscale("log")
        ax.set_title(title)
        ax.grid(True, alpha=0.3)

    # Place legends side by side across the top.
    series_legend = fig.legend(
        series_handles.values(),
        series_handles.keys(),
        loc="upper center",
        ncol=3,
        bbox_to_anchor=(0.42, 0.95),
        title="Method",
    )
    device_legend = fig.legend(
        handles=[
            Line2D([], [], color="k", linestyle="--", marker="s", label="CPU"),
            Line2D([], [], color="k", linestyle="-", marker="o", label="GPU"),
        ],
        labels=["CPU", "GPU"],
        loc="upper center",
        bbox_to_anchor=(0.72, 0.95),
        ncol=1,
        title="Device",
    )
    fig.add_artist(series_legend)
    fig.add_artist(device_legend)
    fig.tight_layout(rect=(0, 0, 1, 0.8))
    fig.savefig(plot_path, bbox_inches="tight")
    print(f"Saved combined plot to {plot_path}")


def _infer_n_frames(csv_path: Path) -> int | None:
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        if reader.fieldnames and "n_frames" in reader.fieldnames:
            for row in reader:
                try:
                    return int(float(row["n_frames"]))
                except (KeyError, ValueError, TypeError):
                    continue
    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot benchmark elapsed times from CSV.")
    parser.add_argument(
        "csv",
        nargs="?",
        default=Path("results/results.csv"),
        type=Path,
        help="Benchmark results CSV (default: results/results.csv)",
    )
    parser.add_argument("--out", type=Path, help="Output plot path (defaults to CSV stem with .png)")
    parser.add_argument(
        "--n-frames",
        type=int,
        help="Number of frames (if absent, inferred from CSV n_frames column; fallback 1000)",
    )
    args = parser.parse_args()

    csv_path = args.csv
    plot_path = args.out if args.out is not None else csv_path.with_suffix(".png")
    plot_path.parent.mkdir(parents=True, exist_ok=True)

    n_frames = args.n_frames if args.n_frames is not None else _infer_n_frames(csv_path) or 1000
    plot_elapsed(csv_path, plot_path, n_frames=n_frames)
