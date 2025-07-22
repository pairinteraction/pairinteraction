# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

# In addition to the `pairinteraction` command installed via the pyproject.toml file,
# this file allows to run the pairinteraction cli via `python -m pairinteraction`.

from pairinteraction.cli import main

if __name__ == "__main__":
    main()
