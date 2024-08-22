#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2024 Sebastian Weber, Henri Menke, and contributors. All rights reserved.
# SPDX-License-Identifier: GPL-3.0-or-later
import runpy


def main():
    runpy.run_module("pairinteraction.gui.app", {}, "__main__")


if __name__ == "__main__":
    main()
