# Pairinteraction - A Rydberg Interaction Calculator

[![Python Wheel](https://github.com/pairinteraction/pairinteraction/actions/workflows/python-wheel.yml/badge.svg)](https://github.com/pairinteraction/pairinteraction/actions/workflows/python-wheel.yml)
[![C++ Coverage - ctest][coverage-cpp-ctest-svg]][coverage-cpp-ctest-link]
[![C++ Coverage - pytest][coverage-cpp-pytest-svg]][coverage-cpp-pytest-link]
[![Python Coverage - pytest][coverage-python-pytest-svg]][coverage-python-pytest-link]
[![PyPI Package][pypi-svg]][pypi-link]
[![arXiv:1612.08053][arXiv-svg]][arXiv-link]
[![License: LGPL v3][license-lgpl-svg]][license-lgpl-link]
[![License: GPL v3][license-gpl-svg]][license-gpl-link]

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

The pairinteraction library without the graphical user interface
is licensed under [LGPL v3][license-lgpl-link].
The graphical user interface as well as the combined work of the pairinteraction
library and the graphical user interface is licensed under [GPL v3][license-gpl-link].


[pypi-svg]: https://img.shields.io/pypi/v/pairinteraction.svg?color=orange
[pypi-link]: https://pypi.org/project/pairinteraction/
[coverage-cpp-ctest-svg]: https://img.shields.io/badge/C%2B%2B_coverage-ctest-blue.svg?style=flat
[coverage-cpp-ctest-link]: https://cuddly-adventure-1w1n2vp.pages.github.io/coverage/cpp-ctest/html/index.html
[coverage-cpp-pytest-svg]: https://img.shields.io/badge/C%2B%2B_coverage-pytest-blue.svg?style=flat
[coverage-cpp-pytest-link]: https://cuddly-adventure-1w1n2vp.pages.github.io/coverage/cpp-pytest/html/index.html
[coverage-python-pytest-svg]: https://img.shields.io/badge/Python_coverage-pytest-blue.svg?style=flat
[coverage-python-pytest-link]: https://cuddly-adventure-1w1n2vp.pages.github.io/coverage/python-pytest/html/index.html
[arXiv-svg]: https://img.shields.io/badge/arXiv-1612.08053-b31b1b.svg?style=flat
[arXiv-link]: https://arxiv.org/abs/1612.08053
[license-lgpl-svg]: https://img.shields.io/badge/License-LGPL_v3-blue.svg?style=flat
[license-gpl-svg]: https://img.shields.io/badge/License-GPLv3-blue.svg?style=flat
[license-lgpl-link]: https://www.gnu.org/licenses/lgpl-3.0.html
[license-gpl-link]: https://www.gnu.org/licenses/gpl-3.0.html
[journal-link]: https://doi.org/10.1088/1361-6455/aa743a
