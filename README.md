# Pairinteraction - Calculating Properties of Rydberg Atoms

[![arXiv:1612.08053][arXiv-svg]][arXiv-link]
[![License: LGPL v3][license-lgpl-svg]][license-lgpl-link]
[![CI Workflow][gh-workflow-svg]][gh-workflow-link]
[![C++ Coverage - ctest][coverage-cpp-ctest-svg]][coverage-cpp-ctest-link]
[![C++ Coverage - pytest][coverage-cpp-pytest-svg]][coverage-cpp-pytest-link]
[![Python Coverage - pytest][coverage-python-pytest-svg]][coverage-python-pytest-link]

[arXiv-svg]: https://img.shields.io/badge/arXiv-1612.08053-b31b1b.svg?style=flat
[arXiv-link]: https://arxiv.org/abs/1612.08053
[license-lgpl-svg]: https://img.shields.io/badge/License-LGPL_v3-blue.svg?style=flat
[license-lgpl-link]: https://www.gnu.org/licenses/lgpl-3.0.html
[pypi-svg]: https://img.shields.io/pypi/v/pairinteraction.svg?style=flat
[pypi-link]: https://pypi.org/project/pairinteraction/
[gh-workflow-svg]: https://github.com/atom-pairinteraction/pairinteraction/actions/workflows/python-wheel.yml/badge.svg
[gh-workflow-link]: https://github.com/atom-pairinteraction/pairinteraction/actions/workflows/python-wheel.yml
[coverage-cpp-ctest-svg]: https://img.shields.io/badge/C%2B%2B_coverage-ctest-blue.svg?style=flat
[coverage-cpp-ctest-link]: https://cuddly-adventure-1w1n2vp.pages.github.io/coverage/cpp-ctest/html/index.html
[coverage-cpp-pytest-svg]: https://img.shields.io/badge/C%2B%2B_coverage-pytest-blue.svg?style=flat
[coverage-cpp-pytest-link]: https://cuddly-adventure-1w1n2vp.pages.github.io/coverage/cpp-pytest/html/index.html
[coverage-python-pytest-svg]: https://img.shields.io/badge/Python_coverage-pytest-blue.svg?style=flat
[coverage-python-pytest-link]: https://cuddly-adventure-1w1n2vp.pages.github.io/coverage/python-pytest/html/index.html


**Note: This is a completely new version of the pairinteraction software that is not backward compatible to versions below v1.0. Breaking changes can occur until the software reaches v2.0.**

The *pairinteraction* software calculates properties of Rydberg atoms. The software consists of a Python library and a graphical user interface for obtaining single-atom properties and calculating pair potentials, making use of a high-performance C++ backend. The software can be installed via pip:

```bash
pip install --only-binary pairinteraction --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ pairinteraction
```

You can use the pairinteraction software [as a Python library][tutorial-link], or you can launch its graphical user interface from the command line:

```bash
pairinteraction gui
```

[tutorial-link]: https://www.pairinteraction.org/pairinteraction/sphinx/html/tutorials/tutorials.html


## Highlights

* For calculating Rydberg pair potentials, the software uses a similar approach as the [old version of pairinteraction](https://github.com/pairinteraction/pairinteraction/tree/v0.9.9), the [Alkali.ne Rydberg Calculator](https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator), and the [rydcalc library](https://github.com/ThompsonLabPrinceton/rydcalc). We optimized the construction and diagonalization of Hamiltonians, typically achieving a **speedup of 5-20x** compared to other implementations.

  ![benchmarking results][benchmark-results-link]

  *Figure: Benchmarking the construction and diagonalization of a Hamiltonian of a pair of Rb 60S atoms for 100 different internuclear distances on an AMD Ryzen 7 5700G CPU using Windows 11. The Hilbert space comprises pair states that differ at most by 4 in n, l and 25GHz in energy. When supported, symmetries where used to reduce the Hilbert space size. See the [benchmarking tool][benchmark-tool].*

* The software uses single-channel quantum defect theory (SQDT) and also **multi-channel quantum defect theory (MQDT)** for the accurate description of atoms.

  The construction of Hamiltonians is accelerated by using pre-calculated matrix elements, which are stored in database tables. These tables are automatically downloaded from GitHub \[[1],[2]\] and cached locally when needed.

* **Electric and magnetic fields in arbitrary directions** can be included in the calculations. Diamagnetism is supported.

[benchmark-results-link]: https://github.com/atom-pairinteraction/pairinteraction/tree/master/data/benchmarking_results/0845d67063_1.4.0-cp312-win_amd-ryzen-7-5700g-with-radeon-graphics_reps4.png
[benchmark-tool]: https://github.com/atom-pairinteraction/pairinteraction/tree/master/tools/benchmarking
[1]: https://github.com/pairinteraction/database-sqdt
[2]: https://github.com/pairinteraction/database-mqdt


## How to Cite

If you use pairinteraction in your research, please cite our tutorial paper:

> Sebastian Weber, Christoph Tresp, Henri Menke, Alban Urvoy, Ofer Firstenberg, Hans Peter Büchler, Sebastian Hofferberth, *Tutorial: Calculation of Rydberg interaction potentials*, [J. Phys. B: At. Mol. Opt. Phys. 50, 133001 (2017)][journal-link], [arXiv:1612.08053][arXiv-link]

Pairinteraction relies on quantum defects provided by the community. Consider citing relevant publications for your atomic species alongside pairinteraction.

<p><details>
<summary><b>Click to expand for quantum defect references</b></summary>

| Element | Model                 | Identifier     | References                                                                                                                                                   |
|---------|-----------------------|----------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
| H       | SQDT                  | `H`            | Schrödinger equation for hydrogen                                                                                                                            |
| Li      | SQDT                  | `Li`           | [10.1017/CBO9780511524530] (1994)<br>[10.1103/PhysRevA.34.2889] (1986)                                                                                       |
| Na      | SQDT                  | `Na`           | [10.1088/0953-4075/30/10/009] (1997)<br>[10.1070/QE1995v025n09ABEH000501] (1995)<br>[10.1103/PhysRevA.45.4720] (1992)                                        |
| K       | SQDT                  | `K`            | [10.1088/0031-8949/27/4/012] (1983)<br>[10.1016/0030-4018(81)90225-X] (1981)                                                                                 |
| Rb      | SQDT                  | `Rb`           | [10.1103/PhysRevA.83.052515] (2011)<br>[10.1103/PhysRevA.74.054502] (2006)<br>[10.1103/PhysRevA.74.062712] (2006)<br>[10.1103/PhysRevA.67.052502] (2003)     |
| Cs      | SQDT                  | `Cs`           | [10.1103/PhysRevA.93.013424] (2016)<br>[10.1103/PhysRevA.35.4650] (1987)<br>[10.1103/PhysRevA.26.2733] (1982)                                                |
| Sr88    | SQDT, singlet sector  | `Sr88_singlet` | [10.1103/PhysRevA.108.022815] (2023)<br>[10.17169/refubium-34581] (2022)                                                                                     |
| Sr88    | SQDT, triplet sector  | `Sr88_triplet` | [10.1016/j.cpc.2020.107814] (2021)                                                                                                                           |
| Sr87    | MQDT                  | `Sr87_mqdt`    | [10.1088/1361-6455/ab4c22] (2019)                                                                                                                            |
| Sr88    | MQDT                  | `Sr88_mqdt`    | [10.1088/1361-6455/ab4c22] (2019)                                                                                                                            |
| Yb171   | MQDT                  | `Yb171_mqdt`   | [10.48550/arXiv.2406.01482] (2024)                                                                                                                           |
| Yb173   | MQDT                  | `Yb173_mqdt`   | MQDT model formulated by us                                                                                                                                  |
| Yb174   | MQDT                  | `Yb174_mqdt`   | [10.48550/arXiv.2406.01482] (2024)                                                                                                                           |

The identifier can be used to specify an atomic species in the pairinteraction software.

</details></p>

[journal-link]: https://doi.org/10.1088/1361-6455/aa743a
[10.1103/PhysRevA.34.2889]: https://doi.org/10.1103/PhysRevA.34.2889
[10.1017/CBO9780511524530]: https://doi.org/10.1017/CBO9780511524530
[10.1103/PhysRevA.45.4720]: https://doi.org/10.1103/PhysRevA.45.4720
[10.1070/QE1995v025n09ABEH000501]: https://doi.org/10.1070/QE1995v025n09ABEH000501
[10.1088/0953-4075/30/10/009]: https://doi.org/10.1088/0953-4075/30/10/009
[10.1088/0031-8949/27/4/012]: https://doi.org/10.1088/0031-8949/27/4/012
[10.1016/0030-4018(81)90225-X]: https://doi.org/10.1016/0030-4018(81)90225-X
[10.1103/PhysRevA.83.052515]: https://doi.org/10.1103/PhysRevA.83.052515
[10.1103/PhysRevA.67.052502]: https://doi.org/10.1103/PhysRevA.67.052502
[10.1103/PhysRevA.74.054502]: https://doi.org/10.1103/PhysRevA.74.054502
[10.1103/PhysRevA.74.062712]: https://doi.org/10.1103/PhysRevA.74.062712
[10.1103/PhysRevA.93.013424]: https://doi.org/10.1103/PhysRevA.93.013424
[10.1103/PhysRevA.26.2733]: https://doi.org/10.1103/PhysRevA.26.2733
[10.1103/PhysRevA.35.4650]: https://doi.org/10.1103/PhysRevA.35.4650
[10.1103/PhysRevA.108.022815]: https://doi.org/10.1103/PhysRevA.108.022815
[10.17169/refubium-34581]: https://doi.org/10.17169/refubium-34581
[10.1016/j.cpc.2020.107814]: https://doi.org/10.1016/j.cpc.2020.107814
[10.1088/1361-6455/ab4c22]: https://doi.org/10.1088/1361-6455/ab4c22
[10.48550/arXiv.2406.01482]: https://doi.org/10.48550/arXiv.2406.01482

## Documentation

**User Guide**

- [Installation] - Whether you prefer an installation via pip or compiling from source, we've got all the information you need.

- [Tutorials] - Covering both the graphical user interface and the Python library,
these tutorials enable you to leverage pairinteraction for your projects.

- [API Reference] - Documentation of classes and functions of pairinteraction's Python library.

- [Publications] - Find out more about the theory behind pairinteraction and how to cite it and the quantum defects you use.

**Contributor Guide**

- [Getting Started as a Contributor] - Discover the many ways you can help improve pairinteraction, from contributing to the repository to providing quantum defects.

- [Overview About pairinteraction's Architecture] - Find out how pairinteraction is structured.

- [Database Format] - Explore the database for storing atomic states and matrix elements.

**Utility Tools [External Links]**

- [jl-mqdt] - Learn how to calculate states and matrix elements using multi-channel quantum defect theory with our tool written in Julia.

- [ryd-numerov] - Learn how to calculate states and matrix elements using single quantum defect theory with our tool written in Python.

[Installation]: https://www.pairinteraction.org/pairinteraction/sphinx/html/installation.html
[Tutorials]: https://www.pairinteraction.org/pairinteraction/sphinx/html/tutorials.html
[API Reference]: https://www.pairinteraction.org/pairinteraction/sphinx/html/reference.html
[Publications]: https://www.pairinteraction.org/pairinteraction/sphinx/html/publications.html
[Getting Started as a Contributor]: https://www.pairinteraction.org/pairinteraction/sphinx/html/getting_started.html
[Overview About pairinteraction's Architecture]: https://www.pairinteraction.org/pairinteraction/sphinx/html/architecture.html
[Database Format]: https://www.pairinteraction.org/pairinteraction/sphinx/html/database.html
[jl-mqdt]: https://github.com/pairinteraction/jl-mqdt/
[ryd-numerov]: https://github.com/pairinteraction/ryd-numerov/


## Contributors

The development of the pairinteraction software has been supported by the [Institute for Theoretical Physics III] of the University of Stuttgart, the Federal Ministry of Education and Research under the Grants [QRydDemo] and [MUNIQC-Atoms], and the company [Atom Computing]. The development of the original version of the software started at the [5th Institute of Physics] of the University of Stuttgart.

The software is maintained by:
* [Sebastian Weber]
* [Johannes Mögerle]

In addition, the following people contributed significantly to the current and/or previous versions of the software:
* [Henri Menke]
* [Frederic Hummel] - Julia package for multi-channel quantum defect theory, matrix elements
* [Eduard J. Braun] - Perturbative calculations, installation instructions for Windows
* [Johannes Block] - Calculation of Rydberg pair potentials near surfaces *(not yet in new version)*
* [Simon Hollerith] - Documentation of the graphical user interface *(not yet in new version)*

We warmly welcome new contributions! Please see our [contributor guide][contributor-link] for more information.

[Institute for Theoretical Physics III]: https://www.itp3.uni-stuttgart.de/
[QRydDemo]: https://www.quantentechnologien.de/forschung/foerderung/quantenprozessoren-und-technologien-fuer-quantencomputer/qryddemo.html
[MUNIQC-Atoms]: https://www.quantentechnologien.de/forschung/foerderung/quantencomputer-demonstrationsaufbauten/muniqc-atoms.html
[Atom Computing]: https://atom-computing.com/
[5th Institute of Physics]: https://www.pi5.uni-stuttgart.de/
[Sebastian Weber]: https://github.com/seweber
[Johannes Mögerle]: https://github.com/johannes-moegerle
[Henri Menke]: https://github.com/hmenke
[Frederic Hummel]: https://github.com/frederic-atom
[Eduard J. Braun]: https://github.com/EduardJBraun
[Johannes Block]: https://github.com/johblock
[Simon Hollerith]: https://github.com/SimonHollerith
[contributor-link]: https://www.pairinteraction.org/pairinteraction/sphinx/html/contribute/ways.html


## License

The pairinteraction software is licensed under [LGPL v3][license-lgpl-link].
