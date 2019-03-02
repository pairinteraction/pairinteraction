# Pairinteraction - A Rydberg Interaction Calculator

[![Travis Build Status][travis-svg]][travis-link]
[![AppVeyor Build Status][appveyor-svg]][appveyor-link]
[![Code Coverage Report][codecov-svg]][codecov-link]
[![arXiv:1612.08053][arXiv-svg]][arXiv-link]
[![License][license-svg]][gpl-link]
   
The *pairinteraction* software calculates properties of Rydberg systems. Visit the official website at https://pairinteraction.github.io/ for documentation and tutorials. Binary builds are available through [GitHub Releases](https://github.com/pairinteraction/pairinteraction/releases).

## Please cite us

> Sebastian Weber, Christoph Tresp, Henri Menke, Alban Urvoy, Ofer Firstenberg, Hans Peter Büchler, Sebastian Hofferberth, *Tutorial: Calculation of Rydberg interaction potentials*, [J. Phys. B: At. Mol. Opt. Phys. 50, 133001 (2017)][journal-link], [arXiv:1612.08053][arXiv-link]

## License

The pairinteraction software consists of the [pairinteraction library](pairinteraction) licensed under the [LGPL v3][lgpl-link], an extension for calculating
radial wave functions using Whittaker functions licensed under the [GPL v3][gpl-link], and the [graphical user interface](pairinteraction_gui) licensed
under the GPL v3. The GPL v3 also applies to the combined work and all provided binary builds.  Note that not all files in this
repository belong to the pairinteraction project. The submodules link to their host repositories. The individual licenses apply.

[travis-svg]: https://travis-ci.org/pairinteraction/pairinteraction.svg?branch=master
[travis-link]: https://travis-ci.org/pairinteraction/pairinteraction
[appveyor-svg]: https://ci.appveyor.com/api/projects/status/t5l4unwt210gq3al/branch/master?svg=true
[appveyor-link]: https://ci.appveyor.com/project/pairinteraction/pairinteraction/branch/master
[codecov-svg]: https://codecov.io/gh/pairinteraction/pairinteraction/branch/master/graph/badge.svg
[codecov-link]: https://codecov.io/gh/pairinteraction/pairinteraction
[arXiv-svg]: doc/sphinx/images/arXiv-badge.svg
[arXiv-link]: https://arxiv.org/abs/1612.08053
[license-svg]: doc/sphinx/images/license-badge.svg
[gpl-link]: https://www.gnu.org/licenses/gpl-3.0.html
[lgpl-link]: https://www.gnu.org/licenses/lgpl-3.0.html
[journal-link]: https://doi.org/10.1088/1361-6455/aa743a
