# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import argparse
import logging
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import duckdb

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s", handlers=[logging.StreamHandler()])

logger = logging.getLogger(__name__)


@dataclass
class TableFile:
    """Dataclass for storing the name, version, and path to a table of states or matrix elements."""

    name: str
    version: int
    path: Path


def get_source_directory() -> Path:
    """Get the source directory of the parquet files."""
    if "PAIRINTERACTION_CACHE_DIR" in os.environ:
        source_dir = Path(os.environ["PAIRINTERACTION_CACHE_DIR"]) / "database"
    elif sys.platform == "win32":
        local_app_data = os.environ.get("LOCALAPPDATA", Path(r"~\AppData\Local").expanduser())
        source_dir = Path(local_app_data) / "pairinteraction" / "database"
    elif sys.platform == "darwin":
        source_dir = Path("~/Library/Caches").expanduser() / "pairinteraction" / "database"
    else:
        xdg_cache = os.environ.get("XDG_CACHE_HOME", "~/.cache")
        source_dir = Path(xdg_cache).expanduser() / "pairinteraction" / "database"
    logger.debug("Source directory determined: %s", source_dir)
    return source_dir


def load_parquet_and_csv_files(path_source: Path) -> dict[str, TableFile]:
    """Load the latest parquet files from the source directory."""
    parquet_and_csv_files: dict[str, TableFile] = {}
    for path in list(path_source.glob("*.parquet")) + list(path_source.glob("*.csv")):
        name, version_str = path.stem.rsplit("_v", 1)
        version = int(version_str)
        if name not in parquet_and_csv_files or version > parquet_and_csv_files[name].version:
            parquet_and_csv_files[name] = TableFile(name, version, path)
            logger.debug("Loaded parquet file: %s", parquet_and_csv_files[name])
    return parquet_and_csv_files


def delete_old_parquet_and_csv_files(path_target: Path, parquet_and_csv_files: dict[str, TableFile]) -> None:
    """Delete old parquet files from the target directory."""
    for path in list(path_target.glob("*.parquet")) + list(path_target.glob("*.csv")):
        name, version_str = path.stem.rsplit("_v", 1)
        version = int(version_str)
        if name in parquet_and_csv_files:
            if version < parquet_and_csv_files[name].version or parquet_and_csv_files[name].path.suffix == ".csv":
                path.unlink()
                logger.info("Deleted csv / old parquet file: %s", path.name)
            elif version > parquet_and_csv_files[name].version:
                raise ValueError(
                    f"Version of the table '{name}' in target directory is higher than in source directory."
                )


def write_parquet_files(
    connection: duckdb.DuckDBPyConnection,
    parquet_and_csv_files: dict[str, TableFile],
    path_target: Path,
    options: dict[str, str],
) -> None:
    """Write parquet files to the target directory. This also updates the path attribute of the parquet files."""
    for parquet_file in parquet_and_csv_files.values():
        parquet_file.path = path_target / f"{parquet_file.name}_v{parquet_file.version}.parquet"
        copy_cmd = (
            f"COPY {parquet_file.name} TO '{parquet_file.path}' "
            f"(FORMAT PARQUET, {', '.join(f'{k} {v}' for k, v in options.items())})"
        )
        connection.execute(copy_cmd)
        logger.debug("Wrote parquet file: %s", parquet_file.path)


def print_metadata(connection: duckdb.DuckDBPyConnection, parquet_files: dict[str, TableFile]) -> None:
    """Print metadata for the given parquet files as a table."""
    headers = ["Name", "Compression", "Rows", "Groups", "Rows/Group", "Size"]

    # Initialize a list to hold each row's data
    table_data = []
    total_size_on_disk = 0.0
    for parquet_file in parquet_files.values():
        if parquet_file.path.suffix == ".csv":
            continue

        metadata = connection.execute(
            "SELECT compression, row_group_num_rows FROM parquet_metadata('?')", parquet_file.path
        ).fetchone()

        file_metadata = connection.execute(
            "SELECT num_rows, num_row_groups FROM parquet_file_metadata('?')", parquet_file.path
        ).fetchone()

        compression = metadata[0] if metadata else "Unknown"
        num_rows = file_metadata[0] / 1e3 if file_metadata else 0  # Convert to thousands
        num_row_groups = file_metadata[1] if file_metadata else 0
        num_rows_per_row_group = (metadata[1] / 1e3) if metadata and metadata[1] else 0  # Convert to thousands
        size_on_disk = parquet_file.path.stat().st_size / 1e6 if parquet_file.path.exists() else 0  # Convert to MB
        total_size_on_disk += size_on_disk

        table_data.append(
            [
                f"{parquet_file.name}_v{parquet_file.version}",
                f"{compression:12s}",
                f"{num_rows:5.0f} k",
                f"{num_row_groups:6d}",
                f"{num_rows_per_row_group:8.0f} k",
                f"{size_on_disk:7.2f} MB",
            ]
        )

    # Add the total size on disk to the table
    headers[0] += f", Total Size: {total_size_on_disk:.2f} MB"

    # Determine column widths based on headers and data
    col_widths = [len(header) for header in headers]
    for row in table_data:
        for idx, item in enumerate(row):
            col_widths[idx] = max(col_widths[idx], len(str(item)))

    # Create a format string for the header and rows
    header_fmt = " | ".join([f"{{:<{w}}}" for w in col_widths])
    separator = "-+-".join(["-" * w for w in col_widths])
    row_fmt = " | ".join([f"{{:<{w}}}" for w in col_widths])

    # Output the table
    logger.info(header_fmt.format(*headers))
    logger.info(separator)
    for row in table_data:
        logger.info(row_fmt.format(*row))


def optimize() -> None:
    """Optimize the parquet database files for performance."""
    parser = argparse.ArgumentParser(
        description="Optimize the parquet database files for performance. By default, the files from the source "
        "directory are overwritten with the optimized files."
    )
    parser.add_argument(
        "--out", type=Path, default=get_source_directory(), help="The output directory for the shrunk parquet files."
    )
    parser.add_argument(
        "--compression",
        type=str,
        choices=["UNCOMPRESSED", "SNAPPY", "ZSTD"],
        default="ZSTD",
        help="The algorithm for compressing the parquet files.",
    )
    args = parser.parse_args()

    path_target = args.out
    compression = args.compression

    path_source = get_source_directory()
    parquet_and_csv_files = load_parquet_and_csv_files(path_source)

    with duckdb.connect(":memory:") as connection:
        # Print metadata of the original parquet files
        logger.info("Metadata of the original parquet files")
        logger.info("")
        print_metadata(connection, parquet_and_csv_files)
        logger.info("")

        # Read parquet files into DuckDB
        for table_file in parquet_and_csv_files.values():
            connection.execute("CREATE TEMP TABLE ? AS SELECT * FROM '?'", table_file.name, table_file.path)
            logger.debug("Loaded table into DuckDB: %s", table_file.name)

        # Delete old parquet files from the target directory to keep it clean
        delete_old_parquet_and_csv_files(path_target, parquet_and_csv_files)

        # Write optimized parquet files
        write_parquet_files(
            connection, parquet_and_csv_files, path_target, {"COMPRESSION": compression, "ROW_GROUP_SIZE": "100_000"}
        )

        # Print metadata of the new parquet files
        logger.info("Metadata of the newly created parquet files")
        logger.info("")
        print_metadata(connection, parquet_and_csv_files)
        logger.info("")


def shrink() -> None:
    """Shrink the parquet database files by filtering and compressing the data."""
    parser = argparse.ArgumentParser(
        description="Shrink the parquet database files by filtering and compressing the data. The default output "
        "directory is the data/database directory of the pairinteraction repository."
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).parent.parent.parent.parent.parent / "data" / "database",
        help="The output directory for the shrunk parquet files.",
    )
    parser.add_argument(
        "--compression",
        type=str,
        choices=["UNCOMPRESSED", "SNAPPY", "ZSTD"],
        default="ZSTD",
        help="The algorithm for compressing the parquet files.",
    )
    parser.add_argument("--min_n", type=int, default=50, help="The minimum value for the quantum number n.")
    parser.add_argument("--max_n", type=int, default=70, help="The maximum value for the quantum number n.")
    parser.add_argument("--max_f", type=int, default=5, help="The maximum value for the quantum number f.")
    args = parser.parse_args()

    path_target = args.out
    compression = args.compression
    min_n = args.min_n
    max_n = args.max_n
    max_f = args.max_f

    path_source = get_source_directory()
    parquet_files = load_parquet_and_csv_files(path_source)

    with duckdb.connect(":memory:") as connection:
        # Print metadata of the original parquet files
        logger.info("Metadata of the original parquet files")
        logger.info("")
        print_metadata(connection, parquet_files)
        logger.info("")

        # Read and filter parquet files
        all_species: list[str] = []

        if "wigner" in parquet_files:
            parquet_file = parquet_files["wigner"]
            connection.execute(
                "CREATE TEMP TABLE ? AS SELECT * FROM '?' WHERE f_initial BETWEEN 0 AND ? AND f_final BETWEEN 0 AND ?",
                parquet_file.name,
                parquet_file.path,
                max_f,
                max_f,
            )
            logger.debug("Filtered 'wigner' table with max_f=%s", max_f)

        for parquet_file in parquet_files.values():
            if parquet_file.name.endswith("_states"):
                species = parquet_file.name.rsplit("_", 1)[0]
                all_species.append(species)
                connection.execute(
                    "CREATE TEMP TABLE ? AS SELECT * FROM '?' WHERE nu BETWEEN ? AND ? AND f BETWEEN 0 AND ?",
                    parquet_file.name,
                    parquet_file.path,
                    min_n,
                    max_n,
                    max_f,
                )
                logger.debug(
                    "Filtered '%s' with nu between %s and %s and f <= %s", parquet_file.name, min_n, max_n, max_f
                )

        for parquet_file in parquet_files.values():
            if not parquet_file.name.endswith("_states") and parquet_file.name != "wigner":
                species = next((s for s in all_species if parquet_file.name.startswith(s)), "")
                if species:
                    state_db_name = f"{species}_states"
                    connection.execute(
                        "CREATE TEMP TABLE ? AS "
                        "WITH s AS (SELECT * FROM '?') "
                        "SELECT s.* FROM s "
                        "JOIN ? AS s1 ON s.id_initial = s1.id "
                        "JOIN ? AS s2 ON s.id_final = s2.id",
                        parquet_file.name,
                        parquet_file.path,
                        state_db_name,
                        state_db_name,
                    )
                    logger.debug("Shrunk table '%s' based on '%s'", parquet_file.name, state_db_name)

        # Remove parquet files for which a recent version exists in the target directory
        for path in list(path_target.glob("*.parquet")) + list(path_target.glob("*.csv")):
            name, version_str = path.stem.rsplit("_v", 1)
            version = int(version_str)
            if name in parquet_files and version > parquet_files[name].version:
                parquet_files.pop(name)

        # Delete old parquet files from the target directory to keep it clean
        delete_old_parquet_and_csv_files(path_target, parquet_files)

        # Write shrunk parquet files with compression
        write_parquet_files(
            connection, parquet_files, path_target, {"COMPRESSION": compression, "ROW_GROUP_SIZE": "100_000"}
        )

        # Print metadata of the new parquet files
        logger.info("Metadata of the newly created parquet files")
        logger.info("")
        print_metadata(connection, parquet_files)
        logger.info("")
