# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import argparse
import shutil
import sys

from colorama import Fore, Style

from pairinteraction import __version__, configure_logging


def main() -> int:
    """Entry point for the pairinteraction CLI."""
    parser = argparse.ArgumentParser(
        description="Pairinteraction CLI: launch the GUI, run tests, or manage the cache",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  pairinteraction --log-level INFO gui\n"
            "  pairinteraction --log-level INFO test\n"
            "  pairinteraction download Rb Cs\n"
            "  pairinteraction purge"
        ),
    )
    parser.add_argument("--version", action="version", version=f"pairinteraction v{__version__}")
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
    download_parser.add_argument("species", nargs="+", help="list of species to download data for")
    download_parser.set_defaults(func=lambda args: download_databases(args.species))

    # Purge command
    purge_parser = subparsers.add_parser("purge", help="delete all cached data")
    purge_parser.set_defaults(func=lambda _args: purge_cache())

    args = parser.parse_args()

    configure_logging(args.log_level)
    if not args.command:
        parser.print_help()
        return 0

    args.func(args)
    return 0


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


def download_databases(species_list: list[str]) -> int:
    """Download the required data files for the specified species."""
    import pairinteraction.real as pi
    from pairinteraction._backend import get_cache_directory

    database_dir = get_cache_directory() / "database"
    database = pi.Database(download_missing=True, use_cache=False, database_dir=database_dir)

    is_wigner_downloaded = False
    exit_code = 0

    for species in species_list:
        print(f"Downloading tables for {species}...")

        try:
            # We make use of the fact that all tables of a species get downloaded
            # automatically when we create a BasisAtom object.
            basis = pi.BasisAtom(species, n=(50, 51), l=(0, 2), database=database)

            # We calculate matrix elements to ensure that the Wigner table is
            # downloaded as well.
            if not is_wigner_downloaded:
                basis.get_matrix_elements(basis, "electric_dipole", 0)
                is_wigner_downloaded = True

            print(Fore.GREEN + "Download successful." + Style.RESET_ALL)
        except Exception as e:
            exit_code = 1
            print(Fore.RED + f"Download failed: {e}" + Style.RESET_ALL)

    return exit_code


def purge_cache() -> int:
    """Delete all cached data."""
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
