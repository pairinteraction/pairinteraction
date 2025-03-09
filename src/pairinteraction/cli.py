import argparse
import logging
import sys

from colorama import Fore, Style, just_fix_windows_console

from pairinteraction import __version__


def main() -> int:
    """Entry point for the pairinteraction CLI."""
    parser = argparse.ArgumentParser(
        description="Pairinteraction CLI: download database tables, launch the GUI, or run tests",
        epilog="Example: pairinteraction --log-level INFO gui",
    )
    parser.add_argument("--version", action="version", version=f"pairinteraction v{__version__}")
    parser.add_argument(
        "--log-level",
        default="WARNING",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="set the logging level (default: WARNING)",
    )
    subparsers = parser.add_subparsers(dest="command", help="available commands", required=True)

    # Download command
    download_parser = subparsers.add_parser("download", help="download database tables for one or more species")
    download_parser.add_argument("species", nargs="+", help="list of species to download data for")
    download_parser.set_defaults(func=lambda args: download_databases(args.species))

    # GUI command
    gui_parser = subparsers.add_parser("gui", help="launch the graphical user interface")
    gui_parser.set_defaults(func=lambda args: start_gui())

    # Test command
    test_parser = subparsers.add_parser("test", help="run module tests")
    test_parser.set_defaults(func=lambda args: run_module_tests())

    args = parser.parse_args()

    setup_logging(args.log_level)

    return args.func(args)


def setup_logging(level_str: str) -> None:
    """Set up logging with a consistent format."""

    class ColoredFormatter(logging.Formatter):
        COLORS = {
            "DEBUG": Fore.BLUE,
            "INFO": Fore.GREEN,
            "WARNING": Fore.YELLOW,
            "ERROR": Fore.RED,
            "CRITICAL": Fore.RED + Style.BRIGHT,
        }

        def format(self, record: logging.LogRecord) -> str:
            record.levelname = f"{self.COLORS.get(record.levelname, '')}{record.levelname}{Style.RESET_ALL}"
            return super().format(record)

    just_fix_windows_console()

    handler = logging.StreamHandler()
    handler.setFormatter(
        ColoredFormatter(
            "%(levelname)s [%(asctime)s.%(msecs)03d] [%(filename)s:%(lineno)d] %(message)s", datefmt="%H:%M:%S"
        )
    )

    logging.basicConfig(level=getattr(logging, level_str.upper()), handlers=[handler])


def download_databases(species_list: list[str]) -> int:
    """Download the required data files for the specified species."""
    import pairinteraction.real as pi

    pi.Database.initialize_global_database(download_missing=True)
    database = pi.Database.get_global_database()

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
                basis.get_matrix_elements(basis, "ELECTRIC_DIPOLE", 0)
                is_wigner_downloaded = True

            print(Fore.GREEN + "Download successful." + Style.RESET_ALL)
        except Exception as e:
            exit_code = 1
            print(Fore.RED + f"Download failed: {e}" + Style.RESET_ALL)

    return exit_code


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
    exit_code = run_module_tests()
    if exit_code:
        print(Fore.RED + "Tests failed." + Style.RESET_ALL)
    else:
        print(Fore.GREEN + "Tests passed." + Style.RESET_ALL)
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
