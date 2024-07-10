# Pairinteraction - A Rydberg Interaction Calculator

[![Linux](https://github.com/pairinteraction/pairinteraction/actions/workflows/linux.yml/badge.svg)](https://github.com/pairinteraction/pairinteraction/actions/workflows/linux.yml)
[![Windows](https://github.com/pairinteraction/pairinteraction/actions/workflows/windows.yml/badge.svg)](https://github.com/pairinteraction/pairinteraction/actions/workflows/windows.yml)
[![macOS](https://github.com/pairinteraction/pairinteraction/actions/workflows/macos.yml/badge.svg)](https://github.com/pairinteraction/pairinteraction/actions/workflows/macos.yml)
[![Coverage Report][codecov-svg]][codecov-link]
[![PyPI Package][pypi-svg]][pypi-link]
[![arXiv:1612.08053][arXiv-svg]][arXiv-link]
[![License][license-svg]][gpl-link]

The *pairinteraction* software calculates properties of Rydberg systems. Visit the official website at https://www.pairinteraction.org/ for documentation and tutorials.

Binary builds are available through [GitHub Releases](https://github.com/pairinteraction/pairinteraction/releases). For using pairinteraction as a Python 3 library, we recommend the installation via pip by calling

```bash
pip install pairinteraction
```

If pairinteraction was installed via pip, the graphical user interface can be started by executing `start_pairinteraction_gui` from the command line.

## Please cite

> Sebastian Weber, Christoph Tresp, Henri Menke, Alban Urvoy, Ofer Firstenberg, Hans Peter Büchler, Sebastian Hofferberth, *Tutorial: Calculation of Rydberg interaction potentials*, [J. Phys. B: At. Mol. Opt. Phys. 50, 133001 (2017)][journal-link], [arXiv:1612.08053][arXiv-link]

The pairinteraction software relies on quantum defects provided by the community.
Please consider citing relevant publications for your atomic species alongside pairinteraction:
*Rb* [add link to copy a bibtex entry], TODO add all species with the name that we use within the software.

## Credits

TODO: Update credits. Maintainers, list of contributors.

## License

The pairinteraction library is licensed under the [LGPL v3][lgpl-link].
The extension for the graphical user interface is licensed under the [GPL v3][gpl-link].
The GPL v3 also applies to the combined work and all provided binary builds.


[pypi-svg]: https://img.shields.io/pypi/v/pairinteraction.svg?color=orange
[pypi-link]: https://pypi.org/project/pairinteraction/
[codecov-svg]: https://img.shields.io/badge/code-coverage-blue.svg?style=flat
[codecov-link]: https://www.pairinteraction.org/pairinteraction/coverage/html/index.html
[arXiv-svg]: https://img.shields.io/badge/arXiv-1612.08053-b31b1b.svg?style=flat
[arXiv-link]: https://arxiv.org/abs/1612.08053
[license-svg]: https://img.shields.io/badge/License-GPLv3-blue.svg?style=flat
[gpl-link]: https://www.gnu.org/licenses/gpl-3.0.html
[lgpl-link]: https://www.gnu.org/licenses/lgpl-3.0.html
[journal-link]: https://doi.org/10.1088/1361-6455/aa743a
