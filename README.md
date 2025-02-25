# Pairinteraction - Calculating Properties of Rydberg Atoms

[![arXiv:1612.08053][arXiv-svg]][arXiv-link]
[![License: LGPL v3][license-lgpl-svg]][license-lgpl-link]

**Note: This is a completely new version of the pairinteraction software that is not backward compatible to versions bellow v1.0. Breaking changes can occur until the software reaches v2.0.**

The *pairinteraction* software calculates properties of Rydberg atoms. The software consists of a Python library and a graphical user interface for obtaining single-atom properties and calculating pair potentials, making use of a high-performance C++ backend. The software can be installed via pip:

```bash
pip install --only-binary pairinteraction --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ pairinteraction
```

The graphical user interface can be started by executing ``start_pairinteraction_gui`` from the command line.

## Highlights

* For calculating Rydberg pair potentials, the software uses a similar approach as the [old version of pairinteraction](https://github.com/pairinteraction/pairinteraction/tree/v0.9.9), the [Alkali.ne Rydberg Calculator](https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator), and the [rydcalc library](https://github.com/ThompsonLabPrinceton/rydcalc). We optimized the construction and diagonalization of Hamiltonians, typically achieving a **speedup of 5-20x** compared to other implementations.

  ![benchmarking results](data/benchmarking_results/0845d67063_1.3.0-cp312-linux_12th-gen-intel(r)-core(tm)-i7-1260p_reps4.png "Benchmarking results")

  *Figure: Benchmarking the construction and diagonalization of a Hamiltonian of a pair of Rb 60S atoms for 100 different internuclear distances on a i7-1260P CPU. The Hilbert space comprises pair states that differ at most by 4 in n, l and 25GHz in energy. When supported, symmetries where used to reduce the Hilbert space size. See the [benchmarking tool](tools/benchmarking).*

* The software uses single-channel quantum defect theory (SQDT) and also **multi-channel quantum defect theory (MQDT)** for the accurate description of atoms.

* **Electric and magnetic fields in arbitrary directions** can be included in the calculations. Diamagnetism is supported.

## Please Cite

> Sebastian Weber, Christoph Tresp, Henri Menke, Alban Urvoy, Ofer Firstenberg, Hans Peter Büchler, Sebastian Hofferberth, *Tutorial: Calculation of Rydberg interaction potentials*, [J. Phys. B: At. Mol. Opt. Phys. 50, 133001 (2017)][journal-link], [arXiv:1612.08053][arXiv-link]

The pairinteraction software relies on quantum defects provided by the community.
Please consider citing relevant publications for your atomic species alongside pairinteraction.

<details>
<summary>Click to expand the table of quantum defect references</summary>

| Species | Model                 | Identifier     | References                                                                                                                                                   |
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

</details>

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

## License

The pairinteraction software is licensed under [LGPL v3][license-lgpl-link].

[arXiv-svg]: https://img.shields.io/badge/arXiv-1612.08053-b31b1b.svg?style=flat
[arXiv-link]: https://arxiv.org/abs/1612.08053
[license-lgpl-svg]: https://img.shields.io/badge/License-LGPL_v3-blue.svg?style=flat
[license-lgpl-link]: https://www.gnu.org/licenses/lgpl-3.0.html
[journal-link]: https://doi.org/10.1088/1361-6455/aa743a
