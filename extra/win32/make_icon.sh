# SPDX-FileCopyrightText: 2024 Sebastian Weber, Henri Menke, and contributors. All rights reserved.
# SPDX-License-Identifier: GPL-3.0-or-later OR LGPL-3.0-or-later

pdflatex icon.tex
convert -density 600 icon.pdf -define icon:auto-resize=256,48,32,16 pairinteraction.ico
