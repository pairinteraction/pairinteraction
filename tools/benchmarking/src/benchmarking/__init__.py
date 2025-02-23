import argparse
import json
import logging
import platform
import sys
from collections.abc import Generator
from contextlib import contextmanager
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from time import perf_counter_ns
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from cpuinfo import get_cpu_info

import pairinteraction
import pairinteraction.backend.complex as pi_complex
import pairinteraction.backend.real as pi_real

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s", handlers=[logging.StreamHandler()])


@dataclass
class BenchmarkResult:
    """Dataclass for storing benchmark results."""

    operation: str
    software: str
    duration_in_sec: float


@contextmanager
def timer() -> Generator[Callable[[], float], None, None]:
    """Timer context manager."""
    start = perf_counter_ns()
    yield lambda: (perf_counter_ns() - start) / 1e9


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
    software_name = f"pairinteraction, {name}"

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


def plot_results(all_results: list[BenchmarkResult], output: str) -> None:
    """Plot the benchmark results."""
    df = pd.DataFrame(all_results)

    sns.set_theme(style="ticks", rc={"axes.spines.right": False, "axes.spines.top": False})

    unique_softwares = df["software"].unique()
    palette = dict(zip(unique_softwares, sns.color_palette("viridis", len(unique_softwares))))

    ax = sns.barplot(
        x="operation", y="duration_in_sec", hue="software", data=df, palette=palette, errorbar="sd", capsize=0.1
    )
    ax.minorticks_off()

    ax.set_xlabel("")
    ax.set_ylabel("Duration (s)")
    ax.legend(title="", ncols=1, loc="upper left")

    for p in ax.patches:
        height = p.get_height()
        if height > 0:
            ax.annotate(
                f"{height:.2f}",
                xy=(p.get_x() + p.get_width() / 2, height),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    plt.tight_layout()
    plt.savefig(output)
    plt.close()

    logging.info(f"Plot saved to '{output}'")


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
    args = parser.parse_args()

    settings = {
        "species": "Rb",
        "n": 60,
        "l": 0,
        "j": 0.5,
        "m": 0.5,
        "delta_n": 3,
        "delta_l": 3,
        "delta_energy_in_GHz": 25,
        "energy_range_in_GHz": 2.5,
        "distances_in_um": np.linspace(1, 10, 100).tolist(),
        "order": 3,
        "download_missing": args.download_missing,
        "diagonalizer": "Lapacke",
    }

    all_results: list[BenchmarkResult] = []

    # Benchmark pairinteraction
    backends = {
        "complex 64": [pi_complex, "float64"],
        "complex 32": [pi_complex, "float32"],
        "real 64": [pi_real, "float64"],
        "real 32": [pi_real, "float32"],
    }
    for _ in range(args.reps):
        for name, [module, float_type] in backends.items():
            logging.info(f"Benchmarking 'pairinteraction, {name}'")
            all_results += benchmark_pairinteraction(module, float_type, name, settings)

    # Generate meaningful output filenames
    hashed = sha256()
    hashed.update(json.dumps(settings, sort_keys=True).encode())
    system = {"Linux": "linux", "Windows": "win", "Darwin": "macosx"}.get(platform.system(), "unknown")
    identifier = f"{pairinteraction.__version__}-cp{sys.version_info.major}{sys.version_info.minor}-{system}"
    cpuname = "-".join(get_cpu_info().get("brand_raw", "unknown").lower().split())

    outdir = Path(__file__).parent.parent.parent.parent.parent / "data" / "benchmarking_results"
    plot_path = outdir / f"{hashed.hexdigest()[:10]}_{identifier}_{cpuname}_reps{args.reps}.png"
    settings_path = outdir / f"{hashed.hexdigest()[:10]}.json"

    # Save plot
    plot_results(all_results, plot_path)

    # Save settings
    with settings_path.open("w", encoding="utf-8") as f:
        json.dump(settings, f, indent=4)
