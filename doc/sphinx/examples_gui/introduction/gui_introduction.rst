work in progress. content of the file will change.

Quantum numbers and atomic species can be specified easily. 

Matrix elements to different principal quantum numbers decrease to to smaller spatial overlaps, the number of different n to account for depends on precision required

Delta l,j, mj which have to be accounted depend on multipol terms due to selection rules for different multipole matrix elements.

The energy bandwidth has to be adjusted according to the principal quantum number (n<-3) because at higher principal quantum number the density of Rydberg states per energy band increases.  At the same time, interactions are more long-range for higher n.

The larger the distances, the less pair states have to be included into diagonalization. For van der Waals potentials, one typically only needs a few pair states.

One can choose whether one wants to only diagonalize in a given symmetry subspace (required for high-precision calculations to save computation time) or calculate all subspaces at the same time (which works for less precise calculations or for the calculation of van der Waals interactions at large distances).

One can choose how the interatomic axis is aligned relative to the z-axis (in the default setting, they are parallel)

External fields not aligned with the interatomic axis break the symmetry of the interaction Hamiltonian and does not allow to do calculations in subspaces.

The plot settings can be adjusted. Colormaps indicate the overlap of the pair potential with the chosen Rydberg pair state. 

The calculated results can be exported and further processed (e.g. using python). The settings can be specified in “Misc”. It can be useful to specify the energy bandwidth where the potentials are exported to save memory.





Vibrational spectroscopy of macrodimer binding potentials:

Macrodimer spectroscopy is currently the most precise way to benchmark Rydberg interactions. (Cite Science paper and recent review?)

Pair interaction can calculate binding potentials at high precision. The config files for a 0+g potential observed blue-detuned from the 35P resonance and a 1_u potential observed blue-detuned from the 36P resonance are included. They can be opened in the gui under “open system settings”

In the config files one can see how the states included in the diagonalization of the interaction Hamiltonian are selected. High-precision calculations require several thousands of pair states. Interestingly, we found the highest precision if the energy bandwidth for the selection of single-atom and two-atom states were similar.

For high precision calculations, symmetries is required to reduce the dimension of the Hamiltonian. For the two exemplary config files we used one gerade and one ungerade potential. Furthermore, pair states are selected according to their angular momentum projection along the interatomic axis.

High-precision calculations in the presence external fields which break the molecular symmetry are very difficult. (see how the dimension of the Hamiltonian increases if applying external fields in the attached config files). For such a case it might be easier to export the pair potentials and its electronic states and calculate the effects of the fields using perturbation theory (if the field is small enough)

After the binding potentials are exported, the vibrational energies are calculated by solving the Schrödinger equation in the interatomic motion assuming adiabatic motion in the exported potentials (assuming the Born-Oppenheimer approximation).



Rydberg interactions have their strongest contribution from dipole-dipole interactions, followed by dipole-quadrupole interactions. Predicting vibrational energies of macrodimers requires to also include higher order multipole terms. In order to approach the experimental accuracy, one needs multipole terms up to 1/R6 (octupole-octupole and other terms scaling similarly). 

The further away the binding potentials from the asymptote, the more pair states are required.

It is just coincidence that both exemplary potentials have another pair potential crossing close to the binding potential minimum




Stark maps.

For Stark maps, only single atom properties have to be specified. In the pair interaction software, this is done by using the same quantum numbers for both atoms.

The number of single-atom states included in the Stark maps can be specified by the single-atom bandwidth. The two-atom energy bandwidth has no meaning in the Stark maps because Stark maps are single-atom properties.


Calculating a Stark map requires to choose different E-fields at the start and the end of the calculation. 

As for the interactions, one can easily specify the quantum numbers of other Rydberg states accounted for.

Stark maps can also be calculated for applied magnetic fields. Here, one chooses the same magnetic field for at the start and at the end

One can choose the relative orientation between the electric and magnetic field. 

One finds that the Stark maps, in particular at small background electric fields such as the ones present in many experiments, depends on the relative orientation between electric and magnetic field. As a consequence, performing Stark spectroscopy or microwave spectroscopy of neighboring Rydberg states can be used to not only measure the electric field strength but also its orientation relative to the applied magnetic bias field.

