import argparse
import json
import logging
import platform
import subprocess
import sys
from collections.abc import Iterable
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Callable

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns
from cpuinfo import get_cpu_info

import pairinteraction.complex as pi_complex
import pairinteraction.real as pi_real
from benchmarking.timer import timer
from pairinteraction import __version__, configure_logging

configure_logging("INFO")


@dataclass
class BenchmarkResult:
    """Dataclass for storing benchmark results."""

    operation: str
    software: str
    duration_in_sec: float


def benchmark_pairinteraction(pi: Callable, float_type: str, name: str, settings: dict) -> list[BenchmarkResult]:
    """Benchmark pairinteraction."""
    species = settings["species"]
    n = settings["n"]
    l = settings["l"]
    j = settings["j"]
    m = settings["m"]
    delta_n = settings["delta_n"]
    delta_l = settings["delta_l"]
    delta_energy = settings["delta_energy_in_GHz"]
    energy_range = settings["energy_range_in_GHz"]
    distances = settings["distances_in_um"]
    order = settings["order"]
    download_missing = settings["download_missing"]
    diagonalizer = settings["diagonalizer"]

    results = []
    software_name = f"{name}"

    with timer() as get_duration:
        db = pi.Database(download_missing=download_missing)
        ket = pi.KetAtom(species, n=n, l=l, j=j, m=m, database=db)
        pair_energy = 2 * ket.get_energy(unit="GHz")
        basis = pi.BasisAtom(species, n=(n - delta_n, n + delta_n), l=(l - delta_l, l + delta_l), database=db)
        system = pi.SystemAtom(basis)
        pair_basis = pi.BasisPair(
            [system, system],
            energy=(pair_energy - delta_energy, pair_energy + delta_energy),
            energy_unit="GHz",
            m=(2 * m, 2 * m),
            product_of_parities="EVEN",
        )
        pair_systems = [
            pi.SystemPair(pair_basis).set_distance(d, unit="micrometer").set_order(order) for d in distances
        ]
        _ = pair_systems[0].get_hamiltonian()
    results.append(BenchmarkResult("Construction", software_name, get_duration()))

    logging.info(f"Number of kets: {pair_systems[0].basis.number_of_states}")

    with timer() as get_duration:
        pi.diagonalize(
            pair_systems,
            diagonalizer=diagonalizer,
            float_type=float_type,
            sort_by_energy=False,
            energy_range=(pair_energy - energy_range, pair_energy + energy_range),
            energy_unit="GHz",
        )

    logging.info(f"Number of states: {pair_systems[0].basis.number_of_states}")

    results.append(BenchmarkResult("Diagonalization", software_name, get_duration()))

    return results


def benchmark_external_script(script_path: Path, name: str, settings: dict) -> list[BenchmarkResult]:
    """Benchmark external script."""
    try:
        proc_result = subprocess.run(  # noqa: S603
            ["uv", "run", "--script", script_path, json.dumps(settings)],  # noqa: S607
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"External script '{script_path}' failed with error:\n{e.stderr}") from e

    output = proc_result.stdout.strip()
    try:
        json_line = next(line for line in reversed(output.splitlines()) if line.startswith("{") and line.endswith("}"))
    except StopIteration as e:
        raise RuntimeError(f"Unable to extract benchmarking results. Output of the subprocess: '{output}'") from e
    data = json.loads(json_line)

    logging.info(f"Number of kets: {data['number_of_kets']}")
    logging.info(f"Number of states: {data['number_of_states']}")

    return [
        BenchmarkResult("Construction", name, data["duration_construction"]),
        BenchmarkResult("Diagonalization", name, data["duration_diagonalization"]),
    ]


def plot_results(all_results: list[BenchmarkResult], output: str) -> None:
    """Plot the benchmark results."""
    df = pd.DataFrame(all_results)

    sns.set_theme(style="ticks", rc={"axes.spines.right": False, "axes.spines.top": False})

    # Define color palette
    unique_softwares = df["software"].unique()
    palette = {
        "Alkali.ne Rydberg Calculator": sns.color_palette("deep")[3],
        "old pairinteraction": sns.color_palette("deep")[1],
    }
    other_softwares = [s for s in unique_softwares if s not in ["Alkali.ne Rydberg Calculator", "old pairinteraction"]]
    palette.update(dict(zip(other_softwares, sns.color_palette("viridis", len(other_softwares)))))

    # Set up a split plot if necessary
    largest = sorted(df["duration_in_sec"])[-1]
    second_largest = sorted(df["duration_in_sec"])[-2]
    is_split = largest > 3 * second_largest

    fig, axes = plt.subplots(
        ncols=1,
        nrows=2 if is_split else 1,
        sharex=True,
        figsize=(8, 2.5),
        gridspec_kw={"height_ratios": [1, 3] if is_split else [1]},
        constrained_layout=True,
    )
    fig.get_layout_engine().set(h_pad=0, hspace=0, rect=[0.10, 0.03, 0.80, 0.92])

    if not isinstance(axes, Iterable):
        axes = [axes]

    # Plot the data
    for ax in axes:
        sns.barplot(
            x="operation",
            y="duration_in_sec",
            hue="software",
            data=df,
            palette=palette,
            errorbar="sd",
            capsize=0.1,
            ax=ax,
        )

    # Annotate the bars with the duration
    for ax in axes:
        ax.set_xlabel("")
        ax.minorticks_off()
        for p in ax.patches:
            height = p.get_height()
            if height == 0:
                continue
            ax.annotate(
                f"{height:.2f}",
                xy=(p.get_x() + p.get_width() / 2, height),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    # Format the axes
    axes[-1].set_ylabel("Duration (s)")
    x_left, x_right = axes[-1].get_xlim()
    axes[-1].set_xlim(x_left, x_right + 2)
    axes[-1].legend(loc="lower right", frameon=False)
    sns.despine(ax=axes[-1])

    if is_split:
        intersection = second_largest * 1.1
        axes[0].set_ylim(largest - intersection * 1 / 3 + intersection * 0.1, largest + intersection * 0.1)
        axes[0].xaxis.set_ticks_position("none")
        axes[0].set_ylabel("")
        axes[-1].set_ylim(0, intersection)

        axes[0].get_legend().remove()

        tick_interval = np.ceil(intersection * 1 / 3 / 2)
        axes[0].yaxis.set_major_locator(ticker.MultipleLocator(tick_interval))
        axes[-1].yaxis.set_major_locator(ticker.MultipleLocator(tick_interval))

        sns.despine(ax=axes[0], bottom=True)

        # Draw split markers
        d = 0.015
        axes[0].plot((-d, +d), (-d * 3, +d * 3), transform=axes[0].transAxes, color="k", clip_on=False)
        axes[-1].plot((-d, +d), (1 - d, 1 + d), transform=axes[-1].transAxes, color="k", clip_on=False)

    # Save the plot
    plt.savefig(output)
    plt.close()


def run() -> None:
    """Run the benchmarking."""
    parser = argparse.ArgumentParser(
        description="Run benchmarks for the pairinteraction software.",
    )
    parser.add_argument(
        "--reps",
        type=int,
        help="Number of repetitions for each benchmark.",
        default=1,
    )
    parser.add_argument(
        "--download-missing",
        action="store_true",
        help="Download missing database files.",
        default=False,
    )
    parser.add_argument(
        "--floats",
        action="store_true",
        help="Run all benchmarks.",
        default=False,
    )
    args = parser.parse_args()

    settings = {
        "species": "Rb",
        "n": 60,
        "l": 0,
        "j": 0.5,
        "m": 0.5,
        "delta_n": 4,
        "delta_l": 4,
        "delta_energy_in_GHz": 25,
        "energy_range_in_GHz": 2.5,
        "distances_in_um": np.linspace(1, 10, 100).tolist(),
        "order": 3,
        "download_missing": args.download_missing,
        "diagonalizer": "lapacke_evr",
    }

    all_results: list[BenchmarkResult] = []

    if args.floats:
        backends = {
            "pairinteraction, real 32": [pi_real, "float32"],
            "pairinteraction, real 64": [pi_real, "float64"],
            "pairinteraction, complex 32": [pi_complex, "float32"],
            "pairinteraction, complex 64": [pi_complex, "float64"],
        }
        external_scripts: dict[str, Path] = {}
    else:
        backends = {
            "pairinteraction": [pi_real, "float32"],
        }
        external_scripts = {
            "old pairinteraction": Path(__file__).parent / "_run_old_pairinteraction.py",
            "Alkali.ne Rydberg Calculator": Path(__file__).parent / "_run_arc.py",
        }

    for _ in range(args.reps):
        # Benchmark pairinteraction backends
        for name, [module, float_type] in backends.items():
            logging.info(f"Benchmarking '{name}'")
            all_results += benchmark_pairinteraction(module, float_type, name, settings)

        # Benchmark ARC and old pairinteraction
        for name, script_path in external_scripts.items():
            logging.info(f"Benchmarking '{name}'")
            all_results.extend(benchmark_external_script(script_path, name, settings))

    # Generate meaningful output filenames
    hashed = sha256()
    hashed.update(json.dumps(settings, sort_keys=True).encode())
    system = {"Linux": "linux", "Windows": "win", "Darwin": "macosx"}.get(platform.system(), "unknown")
    identifier = f"{__version__}-cp{sys.version_info.major}{sys.version_info.minor}-{system}"
    cpuname = "-".join(get_cpu_info().get("brand_raw", "unknown").lower().split())
    outdir = Path.cwd().parent.parent / "data" / "benchmarking_results"
    plot_path = (
        outdir
        / f"{hashed.hexdigest()[:10]}_{identifier}_{cpuname}_reps{args.reps}{'_floats' if args.floats else ''}.png"
    )
    settings_path = outdir / f"{hashed.hexdigest()[:10]}.json"

    # Save plot and settings
    plot_results(all_results, plot_path)
    with settings_path.open("w", encoding="utf-8") as f:
        json.dump(settings, f, indent=4)
    logging.info(f"Plot saved to '{plot_path}'")
    logging.info(f"Settings saved to '{settings_path}'")
