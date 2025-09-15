# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import argparse
import sys
from pathlib import Path
from typing import Callable, cast

from colorama import Fore, Style

from pairinteraction import __version__, configure_logging


def main() -> int:
    """Entry point for the PairInteraction CLI."""
    parser = argparse.ArgumentParser(
        description="PairInteraction CLI: launch the GUI, run tests, or manage the cache",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  pairinteraction --log-level INFO gui\n"
            "  pairinteraction --log-level INFO test\n"
            "  pairinteraction download Rb Cs\n"
            "  pairinteraction download https://github.com/pairinteraction/database-sqdt/releases/download/v1.2/Rb_v1.2.zip\n"
            "  pairinteraction paths\n"
            "  pairinteraction purge"
        ),
    )
    parser.add_argument("--version", action="version", version=f"PairInteraction v{__version__}")
    parser.add_argument(
        "--log-level",
        default="WARNING",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="set the logging level (default: WARNING)",
    )
    subparsers = parser.add_subparsers(dest="command", title="Available Commands")

    # GUI command
    gui_parser = subparsers.add_parser("gui", help="launch the graphical user interface")
    gui_parser.set_defaults(func=lambda _args: start_gui())

    # Test command
    test_parser = subparsers.add_parser("test", help="run module tests")
    test_parser.set_defaults(func=lambda _args: run_module_tests())

    # Download command
    download_parser = subparsers.add_parser("download", help="download database tables for one or more species")
    download_parser.add_argument("species", nargs="+", help="list of species to download data for / list of urls")
    download_parser.set_defaults(func=lambda args: download_databases(args.species))

    # Paths command
    paths_parser = subparsers.add_parser("paths", help="show config and cache directories")
    paths_parser.set_defaults(func=lambda _args: show_paths())

    # Purge command
    purge_parser = subparsers.add_parser("purge", help="delete all cached data")
    purge_parser.set_defaults(func=lambda _args: purge_cache())

    args = parser.parse_args()

    configure_logging(args.log_level)
    if not args.command:
        parser.print_help()
        return 0

    return cast("Callable[[argparse.Namespace], int]", args.func)(args)


def start_gui() -> int:
    """Launch the graphical user interface."""
    from pairinteraction_gui import main as gui_main

    print("Launching the graphical user interface...")
    gui_main()
    return 0


def run_module_tests() -> int:
    """Run the module tests."""
    from pairinteraction import run_module_tests

    print("Running module tests...")
    exit_code = run_module_tests(download_missing=True)
    if exit_code:
        print(Fore.RED + "Tests failed." + Style.RESET_ALL)
    else:
        print(Fore.GREEN + "Tests passed." + Style.RESET_ALL)
    return exit_code


def _download_database_from_url(url: str, tables_dir: Path) -> int:
    import shutil
    import tempfile
    from urllib.request import urlretrieve
    from zipfile import ZipFile

    from packaging.version import Version

    try:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td) / "tables.zip"

            try:
                msg = f"Downloading {url}..."
                print(msg, end="", flush=True)

                def _hook(blocks: int, block_size: int, total_size: int, _msg: str = msg) -> None:
                    if total_size > 0:
                        pct = min(100, int(blocks * block_size * 100 / total_size))
                        print(f"\r{_msg} {pct:3}%", end="", flush=True)

                urlretrieve(url, tmp, reporthook=_hook)  # noqa: S310
            finally:
                print()

            with ZipFile(tmp) as z:
                roots = {Path(n).parts[0] for n in z.namelist() if Path(n).parts}
                root = roots.pop() if len(roots) == 1 else None

                if not root or f"{root}/" not in z.namelist():
                    raise ValueError("The ZIP archive must contain exactly one top-level folder.")  # noqa: TRY301
                species, version_str = root.rsplit("_v", 1)
                Version(version_str)  # validate version

                to_delete = list(tables_dir.glob(f"{species}*"))
                confirmation = input(
                    f"Do you want delete the tables in {', '.join(p.name for p in to_delete)} and "
                    "replace them with the downloaded tables? (y/N): "
                )
                if confirmation.lower() not in ["y", "yes"]:
                    print(Fore.YELLOW + "Aborted replacing tables." + Style.RESET_ALL)
                    return 1
                for p in to_delete:
                    shutil.rmtree(p)

            shutil.unpack_archive(tmp, tables_dir, format="zip")

    except Exception as e:
        print(Fore.RED + f"Failed: {e}" + Style.RESET_ALL)
        return 1

    else:
        print(Fore.GREEN + "Successful." + Style.RESET_ALL)
        return 0


def download_databases(species_list: list[str]) -> int:
    """Download the required data files for the specified species."""
    from urllib.parse import urlparse

    import pairinteraction as pi
    from pairinteraction._backend import get_cache_directory

    database_dir = get_cache_directory() / "database"
    tables_dir = database_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    database = pi.Database(download_missing=True, use_cache=False, database_dir=database_dir)

    is_wigner_downloaded = False
    exit_code = 0

    for species in species_list:
        # If species is a URL, download and unzip to database/tables
        if urlparse(species).scheme in {"http", "https"}:
            print("Check for tables...")
            exit_code |= _download_database_from_url(species, tables_dir)
            continue

        try:
            print(f"Check for tables for {species}...")

            # We make use of the fact that all tables of a species get downloaded
            # automatically when we create a BasisAtom object.
            basis = pi.BasisAtom(species, n=(50, 51), l=(0, 2), database=database)

            # We calculate matrix elements to ensure that the Wigner table is
            # downloaded as well.
            if not is_wigner_downloaded:
                basis.get_matrix_elements(basis, "electric_dipole", 0)
                is_wigner_downloaded = True

            print(Fore.GREEN + "Successful." + Style.RESET_ALL)
        except Exception as e:
            exit_code = 1
            print(Fore.RED + f"Failed: {e}" + Style.RESET_ALL)

    return exit_code


def show_paths() -> int:
    """Show config and cache directories."""
    from pairinteraction._backend import get_cache_directory, get_config_directory

    print("Config directory:", get_config_directory())
    print("Cache directory:", get_cache_directory())
    print("Database directory:", get_cache_directory() / "database/tables")
    return 0


def purge_cache() -> int:
    """Delete all cached data."""
    import shutil

    from pairinteraction._backend import get_cache_directory

    cache_dir = get_cache_directory()

    confirmation = input("Are you sure you want to delete all cached data? (y/N): ")
    if confirmation.lower() not in ["y", "yes"]:
        print(Fore.YELLOW + "Aborted deletion of cache." + Style.RESET_ALL)
        return 0

    print("Deleting cached data...")
    try:
        shutil.rmtree(cache_dir)
    except Exception as e:
        print(Fore.RED + f"Error while deleting cache: {e}" + Style.RESET_ALL)
        return 1
    print(Fore.GREEN + "Cache deleted." + Style.RESET_ALL)
    return 0


if __name__ == "__main__":
    sys.exit(main())
