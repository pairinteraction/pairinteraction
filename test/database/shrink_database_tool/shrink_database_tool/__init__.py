import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict

import duckdb


@dataclass
class ParquetFile:
    name: str
    version: int
    path: Path


def main() -> None:
    """Copy and shrink the database files.

    Copy the parquet files from 'PAIRINTERACTION_CACHE_DIR/database' to
    'test/database' and restrict them to basis states around |60S> so that
    the database is small enough to be included in the repository."""

    # Get the database source directory
    if "PAIRINTERACTION_CACHE_DIR" in os.environ:
        path_source = Path(os.environ["PAIRINTERACTION_CACHE_DIR"]) / "database"
    elif sys.platform == "win32":
        path_cache = Path(os.path.expandvars(r"%LOCALAPPDATA%"))
        if not path_cache.is_absolute:
            path_cache = Path.expanduser(r"~\AppData\Local")
        path_source = path_cache / "pairinteraction" / "database"
    elif sys.platform == "darwin":
        path_source = Path.expanduser(r"~/Library/Caches") / "pairinteraction" / "database"
    else:
        path_cache = Path(os.path.expandvars(r"$XDG_CACHE_HOME"))
        if path_cache.is_absolute:
            path_cache = Path.expanduser(r"~/.cache")
        path_source = path_cache / "pairinteraction" / "database"

    # Get the database target directory
    path_target = Path(__file__).parent.parent.parent

    # Get the latest parquet files in the source directory
    parquet_files: Dict[str, ParquetFile] = {}
    for path in path_source.glob("*.parquet"):
        name, version = path.stem.rsplit("_v", 1)
        parquet_files[name] = ParquetFile(name, int(version), path)

    # Read in the parquet files
    connection = duckdb.connect(":memory:")

    if "wigner" in parquet_files:
        parquet_file = parquet_files["wigner"]
        connection.execute(
            f"CREATE TEMP TABLE {parquet_file.name} "
            f"AS SELECT * FROM '{parquet_file.path}' "
            "WHERE f_initial BETWEEN 0 AND 5 AND f_final BETWEEN 0 AND 5"
        )

    all_species = []
    for parquet_file in parquet_files.values():
        if parquet_file.name.rsplit("_", 1)[-1] == "states":
            all_species.append(parquet_file.name.rsplit("_", 1)[0])
            connection.execute(
                f"CREATE TEMP TABLE {parquet_file.name} "
                f"AS SELECT * FROM '{parquet_file.path}' "
                "WHERE exp_nu BETWEEN 55 AND 65 AND f BETWEEN 0 AND 5"
            )

    for parquet_file in parquet_files.values():
        if parquet_file.name.rsplit("_", 1)[-1] != "states" and parquet_file.name != "wigner":
            species = [s for s in all_species if parquet_file.name.startswith(s)][0]
            state_db_name = f"{species}_states"
            connection.execute(
                f"CREATE TEMP TABLE {parquet_file.name} AS "
                f"WITH s AS (SELECT * FROM '{parquet_file.path}') "
                f"SELECT * FROM s "
                f"JOIN {state_db_name} AS s1 ON s.id_initial = s1.id "
                f"JOIN {state_db_name} AS s2 ON s.id_final = s2.id"
            )

    # Print the number of rows of the parquet files
    for parquet_file in parquet_files.values():
        num_rows = connection.execute(f"SELECT COUNT(*) FROM {parquet_file.name}").fetchone()[0]
        print(f"{parquet_file.name}: {num_rows} rows")

    # Write the parquet files to the target directory
    for parquet_file in parquet_files.values():
        connection.execute(
            f"COPY {parquet_file.name} TO "
            f"'{path_target/parquet_file.name}_v0.parquet' "
            "(FORMAT PARQUET, COMPRESSION ZSTD)"
        )

    # Print the total size of the parquet files
    size = sum(path.stat().st_size for path in path_target.glob("*.parquet"))
    print(f"Size of parquet files: {size * 1e-6:.6f} megabytes")
