import logging
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter_ns

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from cpuinfo import get_cpu_info

import pairinteraction
import pairinteraction.backend.complexdouble as pi_complexdouble
import pairinteraction.backend.complexfloat as pi_complexfloat
import pairinteraction.backend.double as pi_double
import pairinteraction.backend.float as pi_float

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s", handlers=[logging.StreamHandler()])


@dataclass
class BenchmarkResult:
    operation: str
    software: str
    duration: float  # s


@contextmanager
def timer():
    start = perf_counter_ns()
    yield lambda: (perf_counter_ns() - start) / 1e9


def benchmark_pairinteraction(
    pi: callable,
    name: str,
    n: int,
    l: int,
    j: float,
    m: float,
    delta_n: int,
    delta_l: int,
    delta_energy: float,
    distances: np.ndarray,
    order: int,
    download_missing: bool,
) -> list[BenchmarkResult]:
    results = []
    software_name = f"pairinteraction, {name}"

    with timer() as get_duration:
        db = pi.Database(download_missing=download_missing)
        ket = pi.KetAtom("Rb", n=n, l=l, j=j, m=m, database=db)
        pair_energy = 2 * ket.get_energy(unit="GHz")
        basis = pi.BasisAtom("Rb", n=(n - delta_n, n + delta_n), l=(l - delta_l, l + delta_l), database=db)
        system = pi.SystemAtom(basis)
        pair_basis = pi.BasisPair(
            [system, system],
            energy=(pair_energy - delta_energy / 1e9, pair_energy + delta_energy / 1e9),
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
            diagonalizer="Lapacke",
            sort_by_energy=False,
            energy_range=(pair_energy - 2.5, pair_energy + 2.5),
            energy_unit="GHz",
        )

    logging.info(f"Number of states: {pair_systems[0].basis.number_of_states}")

    results.append(BenchmarkResult("Diagonalization", software_name, get_duration()))

    return results


def plot_results(all_results: list[BenchmarkResult], output: str) -> None:
    df = pd.DataFrame(all_results)

    sns.set_theme(style="ticks", rc={"axes.spines.right": False, "axes.spines.top": False})

    unique_softwares = df["software"].unique()
    palette = dict(zip(unique_softwares, sns.color_palette("viridis", len(unique_softwares))))

    ax = sns.barplot(x="operation", y="duration", hue="software", data=df, palette=palette, errorbar="sd", capsize=0.1)
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


def run(download_missing=False, repetitions=4) -> None:
    n, l, j, m = 60, 0, 0.5, 0.5
    delta_n = 3
    delta_l = 3
    delta_energy = 25e9  # Hz
    distances = np.linspace(1, 10, 100)  # um
    order = 3

    all_results: list[BenchmarkResult] = []

    # Benchmark pairinteraction
    backends = {
        "complexdouble": pi_complexdouble,
        "complexfloat": pi_complexfloat,
        "double": pi_double,
        "float": pi_float,
    }
    for name, module in backends.items():
        for _ in range(repetitions):
            all_results += benchmark_pairinteraction(
                module, name, n, l, j, m, delta_n, delta_l, delta_energy, distances, order, download_missing
            )

    # Create a meaningful output filename
    cpuname = "-".join(get_cpu_info().get("brand_raw", "unknown").lower().split())
    output = (
        Path(__file__).parent.parent.parent.parent.parent
        / f"data/benchmarking_results/comparison_{cpuname}_v{pairinteraction.__version__}.png"
    )

    # Plot the results
    plot_results(all_results, output)
