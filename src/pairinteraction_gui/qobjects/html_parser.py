# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

HTML_DICT = {
    "Ex": "E<sub>x</sub>",
    "Ey": "E<sub>y</sub>",
    "Ez": "E<sub>z</sub>",
    "Bx": "B<sub>x</sub>",
    "By": "B<sub>y</sub>",
    "Bz": "B<sub>z</sub>",
}

GREEK = [
    "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa", "lambda", "mu",
    "nu", "xi", "omicron", "pi", "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega",
    "Delta",
]  # fmt: skip


def parse_html(text: str) -> str:
    if text in HTML_DICT:
        return HTML_DICT[text]
    if text in GREEK:
        return f"<span>&{text};</span>"
    if text.startswith("_"):
        return f"<sub>{text[1:]}</sub>"
    if text.startswith("^"):
        return f"<sup>{text[1:]}</sup>"

    return text
