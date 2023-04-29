Introduction to pairinteraction GUI
===================================

Here we introduce how the pairinteraction GUI can be used to calculate **interactions between Rydberg atoms** and **energy shifts in the presence of applied fields**. The underlying physics is discussed in our <tutorial>(LINK to your paper).

Rydberg pair interactions
-------------------------

In the image below you can see the window of the user interface. Here, it has been used to calculate van der Waals interactions between two s-state Rydberg atoms in rubidium, at principal quantum numbers $n=70$. 

<inlcude Rydberg_SS_interactions.png>

The <config file>(LINK: settings_ss_Rb.sconf) where all settings for this calculations are stored can be imported via "File -> importing system settings" in the upper left. The different variables which can be chosen by the user are explained in the following.

- The **quantum numbers** $n,l,j,m_J$ of both interacting Rydberg atoms. For van der Waals Rydberg interactions between two identical Rydberg atoms, one uses the same quantum numbers for both atoms. For Rydberg atom pairs interacting via resonant dipole-dipole interactions, one chooses states with $Δ_l = 1$. 
- The potentials are calculated by diagonalizing the Rydberg interaction Hamiltonian at various interatomic distances. The **resolution** can be specified in by number of steps. 
- The distance between the atoms at the **start and the end of the distance interval** where interactions should be calculated. Also the orientation of the atom pairs relativ to the z-axis is tunable (per default, the atoms are parallel to the z-axis). 
- **External fields** $E$ and $B$. For pair interaction calculations where the distance is varied, one typically chooses identical field strength from start to end.
- In the calculation, the interaction Hamiltonian is in **multipole terms**. At large distances, it is usually good enough to use the lowest order multipole expansion where only dipole-dipole interactions are accounted for ($1/R^3$). At closer distances, also higher order terms become important. 
- **Number of states used for the diagonalization**. One can select via energy bands for single-atom states and their quantum numbers as well as via energy bands in the pair states. Because the decreasing spatial overlap for large differences in $n$, only nearby Rydberg states have to be included. The single-atom quantum numbers can be selected according to the selection rules from the multipole-matrix elements (i.e., for dipole-dipole interactions, only states with $Δ_l = 1$ have finite coupling elements). 
- The **symmetries of the interaction Hamiltonian** can be used to restrict the basis to the relevant subspaces and therefore shorten the computation time. This only works in the absense of external fields because fields (if they are not pointing along the interatomic axis) break the symmetry of the interaction Hamiltonian.

The results are plotted in the gui, **plot settings** can be manually changed. The colormaps in the plots indicate the state overlap of the pair potentials with the chosen Rydberg state. The calculations can be exported using the save button. Under **Misc**, the settings for the export can be specified. In order to save memory, it makes sense to confine the exported states to the energy band where one is interested in. The calculated results can then be further processed (e.g. using python). Under configuration, the **cache directory** can be specified.



Vibrational spectroscopy of macrodimer binding potentials:
---------------------------------------------------------

Macrodimers <cite https://pubs.acs.org/doi/10.1021/acs.jpca.2c08454> are micrometer-seized diatomic molecules consisting of Rydberg pairs bound by their interactions. Macrodimer spectroscopy is currently the most precise way to benchmark Rydberg interactions because it provides narrow spectroscopic signatures and probes the interactions at short distances. For rubidium, the binding potentials can be calculated at high accuracy using the pair interaction software. For more complicated atoms such as Sr or Yb, macrodimer spectroscopy may provide valuable insights into the details of their Rydberg interactions.

Here, we include the config files for:

- $0^{+}_g$<LINK to settings_0gp.sconf CITE https://www.science.org/doi/10.1126/science.aaw4150> potential observed blue-detuned from the $35P_{1/2}$ resonance.
- $1_u$<LINK to settings_1u.sconf CITE https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.128.113602>  potential observed blue-detuned from the $36P_{1/2}$ resonance. 

In contrast to interactions at large distances which are usually described by van der Waals interactions, these high-precision calculations require several thousands of pair states. The larger the interaction shift from the non-interacting pair state energies, the larger the required basis. In order to reduce the computational effort, the basis should be restricted to the symmetry subspace of interest. In the two attached config files, we used one gerade and one ungerade potential. Furthermore, pair states are selected according to their angular momentum projection along the interatomic axis. We found the highest precision if the energy bandwidth for the selection of single-atom and two-atom states were similar.

High-precision calculations in the presence external fields which break the molecular symmetry are very difficult (see how the dimension of the Hamiltonian increases if applying external fields in the attached config files). For finite fields, if the field strength is small enough, it might be easier to export the pair potentials and its electronic states and calculate the effects of the fields using perturbation theory. 

The vibrational energies can be calculated by solving the Schrödinger equation of the interatomic motion in the calculated binding potentials. Here, in many cases (but not all), the Born-Oppenheimer approximation holds. 

Multipole interactions
^^^^^^^^^^^^^^^^^^^^^^

Rydberg interactions have their strongest contribution from dipole-dipole interactions $1/R^3$, followed by dipole-quadrupole interactions $1/R^4$. Calculating vibrational energies of macrodimers requires to also include higher order multipole terms. Approaching the experimental accuracy requires one to include multipole terms up to $1/R^6$ (octupole-octupole and other terms with the same scaling). The image below shows how the calcations approach the observed vibrational resonance when the number of multipole terms is increased. 

<Include Multipole_terms.png CITE https://edoc.ub.uni-muenchen.de/30114/>


Stark maps
-----------

- The user interface can also be very helpful to calculate energy shifts of Rydberg states in the presence of applied electric and magnetic fields $E$ and $B$. Here, only single atom properties have to be specified. In the pair interaction software, this is done by using the same quantum numbers for both atoms. 

- The number of single-atom states included in the Stark maps can be specified by the single-atom bandwidth. Again, the quantum numbers of the Rydberg states used for the calculation can be specified. The two-atom energy bandwidth has no meaning in the Stark maps because Stark maps are single-atom properties.

- Also the interatomic distance range varied in the calculation of Rydberg interactions is not included in the calculation of Stark maps. Instead, the E-field is varied at the start and the end of the calculation. 

<Include Stark_maps.png>

- Stark maps can also be calculated in the presence of magnetic fields $B$ (see attached config file LINK to settings_Stark_Maps.sconf), also the relative orientation between $E$ and $B$ can be specified in the GUI. In the experiment, this can be useful to obtain the background electric field and its orientation from spectroscopy between different Rydberg states since the splitting depends sligthly on the orientation between both fields. Stark maps can also be calculated at high fields. 

