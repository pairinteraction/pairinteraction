# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import argparse
import re
import sys
from pathlib import Path
from typing import TYPE_CHECKING, cast

from colorama import Fore, Style

from pairinteraction import __version__, configure_logging

if TYPE_CHECKING:
    from collections.abc import Callable


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """Show default arguments and keep manual line breaks in the epilog."""


def main() -> int:
    """Entry point for the PairInteraction CLI."""
    parser = argparse.ArgumentParser(
        description=("PairInteraction CLI\n\nRun 'pairinteraction' without a command to launch the GUI."),
        formatter_class=HelpFormatter,
        epilog=(
            "Examples:\n"
            "  pairinteraction\n"
            "  pairinteraction --log-level INFO\n"
            "  pairinteraction --log-level INFO test\n"
            "  pairinteraction database list\n"
            "  pairinteraction database download Rb Cs\n"
            "  pairinteraction database download Rb Cs --version 2.0\n"
            "  pairinteraction database download https://github.com/pairinteraction/database-sqdt/releases/download/v2.0/Rb_v2.0.zip\n"
            "  pairinteraction database remove\n"
            "  pairinteraction config reset-gui\n"
            "  pairinteraction config paths\n"
            "\n"
            "Command-specific help:\n"
            "  pairinteraction test --help\n"
            "  pairinteraction database --help\n"
            "  pairinteraction config --help"
        ),
    )
    parser.add_argument("--version", action="version", version=f"PairInteraction v{__version__}")
    parser.add_argument(
        "--reload",
        action="store_true",
        help="launch the GUI with automatic theme reload during development",
    )
    parser.add_argument(
        "--log-level",
        default="WARNING",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="set the logging level",
    )

    # Launch GUI (default action)
    parser.set_defaults(func=lambda args: start_gui(reload=args.reload))

    subparsers = parser.add_subparsers(dest="command", title="commands", metavar="{test,database,config}")

    # Removed launch command
    gui_parser = subparsers.add_parser(
        "gui",
        add_help=False,
    )
    gui_parser.set_defaults(
        func=lambda _args: parser.error(
            "The 'gui' subcommand no longer exists. To launch the GUI, run 'pairinteraction' without a command."
        )
    )

    # Test command
    test_parser = subparsers.add_parser(
        "test",
        formatter_class=HelpFormatter,
        help="run tests",
    )
    test_parser.set_defaults(func=lambda _args: run_unit_tests())

    # Database command group
    database_parser = subparsers.add_parser(
        "database",
        formatter_class=HelpFormatter,
        help="manage and inspect the database",
    )
    database_subparsers = database_parser.add_subparsers(dest="database_command", title="database commands")

    # Database list command
    db_list_parser = database_subparsers.add_parser(
        "list",
        formatter_class=HelpFormatter,
        help="list local and remote database table versions",
    )
    db_list_parser.set_defaults(func=lambda _args: list_databases())

    # Database download command
    db_download_parser = database_subparsers.add_parser(
        "download",
        formatter_class=HelpFormatter,
        help="download database tables for one or more species",
    )
    db_download_parser.add_argument("species", nargs="+", help="list of species to download data for / list of URLs")
    db_download_parser.add_argument(
        "--version",
        metavar="X.Y",
        default=None,
        help="version of the database tables to download, e.g. '2.0' or '2' for the newest tables of major version 2"
        " (by default, the newest compatible tables are downloaded); cannot be combined with URLs",
    )
    db_download_parser.set_defaults(func=lambda args: download_databases(args.species, args.version))

    # Database remove command
    db_remove_parser = database_subparsers.add_parser(
        "remove",
        formatter_class=HelpFormatter,
        help="delete the cached database directory",
    )
    db_remove_parser.set_defaults(func=lambda _args: remove_database_cache())

    database_parser.set_defaults(func=lambda _args: print_help(database_parser))

    # Config command group
    config_parser = subparsers.add_parser(
        "config",
        formatter_class=HelpFormatter,
        help="manage GUI settings and inspect paths",
    )
    config_subparsers = config_parser.add_subparsers(dest="config_command", title="config commands")

    # Config reset gui command
    config_reset_gui_parser = config_subparsers.add_parser(
        "reset-gui",
        formatter_class=HelpFormatter,
        help="delete GUI settings file to restore defaults",
    )
    config_reset_gui_parser.set_defaults(func=lambda _args: reset_gui_settings())

    # Config list paths command
    config_list_paths_parser = config_subparsers.add_parser(
        "paths",
        formatter_class=HelpFormatter,
        help="show config and cache directories",
    )
    config_list_paths_parser.set_defaults(func=lambda _args: show_paths())

    config_parser.set_defaults(func=lambda _args: print_help(config_parser))

    args = parser.parse_args()

    if args.command is not None and args.reload:
        parser.error("--reload can only be used when launching the GUI")

    configure_logging(args.log_level)

    return cast("Callable[[argparse.Namespace], int]", args.func)(args)


def print_help(parser: argparse.ArgumentParser) -> int:
    """Print help."""
    parser.print_help()
    return 0


def start_gui(*, reload: bool = False) -> int:
    """Launch the GUI."""
    from pairinteraction_gui import main as gui_main

    print("Launching the GUI...")
    gui_main(enable_theme_hot_reload=reload)
    return 0


def reset_gui_settings() -> int:
    """Delete the GUI settings file to restore default values."""
    from pairinteraction._backend import get_cache_directory

    settings_file = get_cache_directory() / "gui_settings.ini"

    if not settings_file.exists():
        print("No GUI settings file found. Nothing to delete.")
        return 0

    confirmation = input(f"Are you sure you want to delete the GUI settings file {settings_file}? (y/N): ")
    if confirmation.lower() not in ["y", "yes"]:
        print(Fore.YELLOW + "Aborted deletion of GUI settings." + Style.RESET_ALL)
        return 0

    try:
        settings_file.unlink()
    except Exception as e:
        print(Fore.RED + f"Error while deleting GUI settings file: {e}" + Style.RESET_ALL)
        return 1

    print(Fore.GREEN + "GUI settings deleted. Default values will be used on next launch." + Style.RESET_ALL)
    return 0


def run_unit_tests() -> int:
    """Run the C++ module unit tests."""
    from pairinteraction import run_unit_tests

    print("Running the C++ module unit tests...")
    exit_code = run_unit_tests(download_missing=True)
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

    from pairinteraction._backend import COMPATIBLE_DATABASE_VERSION_MAJOR

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
                version = Version(version_str)  # validate version
                compatible_major = COMPATIBLE_DATABASE_VERSION_MAJOR
                if version.major != compatible_major:
                    raise ValueError(  # noqa: TRY301
                        f"The tables in '{root}' use the database format v{version.major}, but this version of "
                        f"PairInteraction requires v{compatible_major}. They would be ignored after downloading, "
                        "so please download tables of a matching version instead."
                    )

                to_delete = list(tables_dir.glob(f"{species}_v*"))
                if to_delete:
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


_TABLE_ASSET_REGEX = re.compile(r"^(\w+)_v(\d+)\.(\d+)\.zip$")


def _parse_requested_version(version: str) -> tuple[int, int | None]:
    """Parse a requested table version of the form 'X.Y' or 'X' (an optional 'v' prefix is allowed)."""
    from pairinteraction._backend import COMPATIBLE_DATABASE_VERSION_MAJOR

    match = re.fullmatch(r"v?(\d+)(?:\.(\d+))?", version)
    if match is None:
        raise ValueError(f"Invalid version '{version}'. Expected a version like '2.0' or '2'.")

    major = int(match[1])
    if major != COMPATIBLE_DATABASE_VERSION_MAJOR:
        raise ValueError(
            f"The requested tables use the database format v{major}, but this version of PairInteraction "
            f"requires v{COMPATIBLE_DATABASE_VERSION_MAJOR}. They would be ignored after downloading, "
            "so please request tables of a matching version instead."
        )

    return major, int(match[2]) if match[2] is not None else None


def _get_releases_endpoints() -> list[str]:
    """Return the URLs listing the releases of the configured database repositories."""
    import json

    from pairinteraction._backend import get_config_directory

    host = "https://api.github.com"
    repo_paths = [
        "/repos/pairinteraction/database-sqdt/releases/latest",
        "/repos/pairinteraction/database-mqdt/releases/latest",
    ]

    config_file = Path(get_config_directory()) / "database.json"
    if config_file.exists():
        try:
            doc = json.loads(config_file.read_text())
            host = doc["database_repo_host"]
            repo_paths = doc["database_repo_paths"]
        except Exception as e:
            print(Fore.YELLOW + f"Ignoring {config_file}: {e}" + Style.RESET_ALL)

    # Turn the endpoints of single releases (e.g. '.../releases/latest') into endpoints listing all releases
    endpoints = [host.rstrip("/") + path.partition("/releases")[0] + "/releases" for path in repo_paths]
    return list(dict.fromkeys(endpoints))


def _fetch_table_assets(names: list[str]) -> dict[str, dict[tuple[int, int], str]]:
    """Look up the tables that are available in the configured database repositories.

    Args:
        names: The names of the tables of interest (i.e., species identifiers or 'misc').

    Returns:
        A dictionary mapping each table name to the download URLs of its available versions.

    """
    import json
    import os
    from urllib.request import Request, urlopen

    headers = {"Accept": "application/vnd.github+json", "User-Agent": f"pairinteraction/{__version__}"}
    token = os.environ.get("GITHUB_TOKEN")
    if token:
        headers["Authorization"] = f"Bearer {token}"

    assets: dict[str, dict[tuple[int, int], str]] = {name: {} for name in names}
    for endpoint in _get_releases_endpoints():
        request = Request(f"{endpoint}?per_page=100", headers=headers)  # noqa: S310
        with urlopen(request, timeout=30) as response:  # noqa: S310
            releases = json.load(response)

        for asset in (asset for release in releases for asset in release.get("assets", [])):
            match = _TABLE_ASSET_REGEX.match(asset["name"])
            if match is not None and match[1] in assets:
                assets[match[1]][int(match[2]), int(match[3])] = asset["browser_download_url"]

    return assets


def _download_databases_of_version(species_list: list[str], version: str) -> int:
    """Download the tables of the specified version for the specified species."""
    from pairinteraction._backend import get_cache_directory

    try:
        major, minor = _parse_requested_version(version)
    except ValueError as e:
        print(Fore.RED + str(e) + Style.RESET_ALL)
        return 1

    tables_dir = get_cache_directory() / "database" / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    # The 'misc' tables (e.g. the table of Wigner 3j symbols) are needed for calculating matrix elements
    names = list(dict.fromkeys([*species_list, "misc"]))

    print(f"Looking up tables of version {version}...")
    try:
        assets = _fetch_table_assets(names)
    except Exception as e:
        print(Fore.RED + f"Failed to look up the available tables: {e}" + Style.RESET_ALL)
        return 1

    exit_code = 0
    for name, urls in assets.items():
        print(f"Check for tables for {name}...")

        matching = [v for v in urls if v[0] == major and (minor is None or v[1] == minor)]
        if matching:
            exit_code |= _download_database_from_url(urls[max(matching)], tables_dir)
            continue

        available = ", ".join(f"{v[0]}.{v[1]}" for v in sorted(urls))
        message = f"Failed: no tables of version {version} available for {name}."
        message += f" Available versions: {available}." if available else ""
        is_fatal = name != "misc"  # missing misc tables are downloaded automatically when they are needed
        print((Fore.RED if is_fatal else Fore.YELLOW) + message + Style.RESET_ALL)
        exit_code |= int(is_fatal)

    return exit_code


def download_databases(species_list: list[str], version: str | None = None) -> int:
    """Download the required data files for the specified species."""
    from urllib.parse import urlparse

    import pairinteraction as pi
    from pairinteraction._backend import get_cache_directory

    if version is not None:
        if any(urlparse(species).scheme in {"http", "https"} for species in species_list):
            print(Fore.RED + "The --version option cannot be combined with URLs." + Style.RESET_ALL)
            return 1
        return _download_databases_of_version(species_list, version)

    database_dir = get_cache_directory() / "database"
    tables_dir = database_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    database: pi.Database | None = None  # created lazily so that URLs can be used without internet access
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

            if database is None:
                database = pi.Database(download_missing=True, use_cache=False, database_dir=database_dir)

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


def list_databases() -> int:
    """Print a table of local and remote database table versions."""
    from pairinteraction.database import print_database_info

    print_database_info()
    return 0


def show_paths() -> int:
    """Show config and cache directories."""
    from pairinteraction._backend import get_cache_directory, get_config_directory

    print("Config directory:", get_config_directory())
    print("Cache directory:", get_cache_directory())
    print("Database directory:", get_cache_directory() / "database/tables")
    return 0


def remove_database_cache() -> int:
    """Delete the cached database directory."""
    import shutil

    from pairinteraction._backend import get_cache_directory

    database_dir = get_cache_directory() / "database"

    confirmation = input(f"Are you sure you want to delete all downloaded database tables in {database_dir}? (y/N): ")
    if confirmation.lower() not in ["y", "yes"]:
        print(Fore.YELLOW + "Aborted deletion of database directory." + Style.RESET_ALL)
        return 0

    print(f"Deleting cached database directory {database_dir}...")
    try:
        shutil.rmtree(database_dir)
    except Exception as e:
        print(Fore.RED + f"Error while deleting database directory: {e}" + Style.RESET_ALL)
        return 1

    print(Fore.GREEN + "Database directory deleted." + Style.RESET_ALL)
    return 0


if __name__ == "__main__":
    sys.exit(main())
