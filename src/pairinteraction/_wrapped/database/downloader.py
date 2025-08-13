# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import json
import logging
import os
import urllib.request
import zipfile
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Any, Literal, Union

from tqdm import tqdm

logger = logging.getLogger(__name__)


@dataclass
class ReleaseAsset:
    """Information about a release asset."""

    id: int
    full_name: str
    download_url: str

    @property
    def name(self) -> str:
        """Extract the asset name from the full name."""
        return self.full_name.replace(".zip", "")

    @property
    def species(self) -> str:
        """Extract species from the asset name."""
        parts = self.name.split("_")
        if len(parts) < 2:
            raise ValueError(f"Invalid asset name format: {self.name}")
        return "_".join(parts[:-1])

    @property
    def version(self) -> tuple[int, int]:
        """Extract version from the asset name."""
        parts = self.name.split("_")
        if len(parts) < 2:
            raise ValueError(f"Invalid asset name format: {self.name}")
        return parse_version(parts[-1])


class GitHubDownloader:
    """Downloads assets from GitHub releases."""

    def __init__(self, repo: str) -> None:
        """Initialize downloader for a specific repository.

        Args:
            repo: Repository in format 'owner/repo'

        """
        self.repo = repo
        self.base_url = f"https://api.github.com/repos/{repo}"

    @cached_property
    def release_data(self) -> dict[str, Any]:
        """Get information about the latest release."""
        url = f"{self.base_url}/releases/latest"

        headers = {"Accept": "application/vnd.github+json"}
        self._add_common_headers(headers)

        request = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(request, timeout=10) as response:
            return json.loads(response.read())  # type: ignore[no-any-return]

    @cached_property
    def assets(self) -> list[ReleaseAsset]:
        """List all assets for the latest release."""
        assets_data = self.release_data.get("assets", [])
        return [
            ReleaseAsset(id=asset["id"], full_name=asset["name"], download_url=asset["url"]) for asset in assets_data
        ]

    def download_and_extract(
        self, species: str, database_dir: Union[str, os.PathLike[str]] = "."
    ) -> Literal["succeeded", "skipped"]:
        """Download the latest release asset for a specific species.

        Args:
            species: Species to download
            database_dir: Directory to save the downloaded file. Defaults to current directory.

        """
        database_dir = Path(database_dir)

        if not database_dir.exists():
            logger.info("Creating database directory %s", database_dir)
            database_dir.mkdir(parents=True, exist_ok=True)

        # check if local folder already exists
        local_folders = [f for f in database_dir.iterdir() if f.name.startswith(f"{species}_v") and f.is_dir()]
        local_versions = [parse_version(f.name.split("_")[-1]) for f in local_folders]
        local_version = max(local_versions, default=None)

        asset = next((a for a in self.assets if a.species == species), None)
        if asset is None:
            raise ValueError(f"No asset found for species '{species}'")

        if local_version is None:
            pass
        elif asset.version < local_version:
            logger.info(
                "Local version %s is newer than remote version %s. Skipping download.", local_version, asset.version
            )
            return "skipped"
        elif asset.version == local_version:
            logger.info(
                "Local version %s is the same as remote version %s. Skipping download.", local_version, asset.version
            )
            return "skipped"

        output_path = database_dir / asset.full_name
        if output_path.suffix != ".zip":
            raise ValueError(f"Expected a .zip file, got {output_path.suffix}")
        if output_path.exists():
            logger.debug("Removing existing file %s", output_path)
            output_path.unlink(missing_ok=True)

        logger.info("Downloading asset for species '%s' version '%s' to %s", species, asset.version, output_path)
        self._download_asset(asset.id, output_path)

        logger.info("Extracting %s to %s", asset.full_name, database_dir)
        with zipfile.ZipFile(output_path, "r") as zip_ref:
            zip_ref.extractall(database_dir)

        logger.debug("Removing downloaded zip file %s", output_path)
        output_path.unlink(missing_ok=True)
        return "succeeded"

    def _download_asset(self, asset_id: int, output_path: Union[str, os.PathLike[str]]) -> None:
        """Download a specific release asset by ID.

        Args:
            asset_id: The asset ID to download
            output_path: Where to save the downloaded file

        """
        url = f"{self.base_url}/releases/assets/{asset_id}"
        headers = {"Accept": "application/octet-stream"}
        self._add_common_headers(headers)

        request = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(request, timeout=10) as response:
            total_size = int(response.headers.get("Content-Length", 0))
            chunk_size = 1024 * 1024  # 1MB chunks

            with Path(output_path).open("wb") as f, tqdm(total=total_size, unit="iB", unit_scale=True) as pbar:
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    size = f.write(chunk)
                    pbar.update(size)

    def _add_common_headers(self, headers: dict[str, str]) -> None:
        """Add common headers to the request.

        This includes the GitHub token if available, as well as the Api-Version.
        """
        headers["X-GitHub-Api-Version"] = "2022-11-28"

        github_token = os.getenv("GITHUB_TOKEN")
        if github_token:
            headers["Authorization"] = f"Bearer {github_token}"


def parse_version(version_str: str) -> tuple[int, int]:
    """Parse a version string into a tuple of (int, int)."""
    if not version_str.startswith("v"):
        raise ValueError(f"Version string must start with 'v': {version_str}")
    parts = version_str[1:].split(".")
    if len(parts) != 2:
        raise ValueError(f"Version string must be in format 'vX.Y': {version_str}")
    return int(parts[0]), int(parts[1])
