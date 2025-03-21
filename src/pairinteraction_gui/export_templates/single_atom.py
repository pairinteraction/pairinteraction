import matplotlib.pyplot as plt
import numpy as np

import pairinteraction.real as pi

# # # # # State of Interest # # # # #
ket = pi.KetAtom(SPECIES, **QUANTUM_NUMBERS)
print(f"State of interest: {ket}")

# # # # # Basis # # # # #
basis = pi.BasisAtom(ket.species, **QUANTUM_NUMBERS_RESTRICTIONS)
print(str(basis))
print(f" ⇒ Basis consists of {basis.number_of_kets} kets")

# # # # # System # # # # #
steps = STEPS
Ex = np.linspace(EX_MIN, EX_MAX, steps)
Ey = np.linspace(EY_MIN, EY_MAX, steps)
Ez = np.linspace(EZ_MIN, EZ_MAX, steps)
Bx = np.linspace(BX_MIN, BX_MAX, steps)
By = np.linspace(BY_MIN, BY_MAX, steps)
Bz = np.linspace(BZ_MIN, BZ_MAX, steps)

systems = []
for i in range(steps):
    system = pi.SystemAtom(basis)
    system.set_electric_field([Ex[i], Ey[i], Ez[i]], unit="V/cm")
    system.set_magnetic_field([Bx[i], By[i], Bz[i]], unit="G")
    systems.append(system)

pi.diagonalize(systems, **DIAGONALIZE_KWARGS)

energies = [system.get_eigenvalues(unit="GHz") - ket.get_energy("GHz") for system in systems]
overlaps = [system.get_eigenbasis().get_overlaps(ket) for system in systems]


# # # # # Plot # # # # #
fig, ax = plt.subplots()

x_values = X_VALUES
ax.set_xlabel(XLABEL)
ax.set_ylabel("Energy [GHz]")


try:
    ax.plot(x_values, np.array(energies), c="0.9", lw=0.25, zorder=-10)
except ValueError:  # inhomogeneous shape -> no simple line plot possible
    for x, es in zip(x_values, energies):
        ax.plot([x] * len(es), es, c="0.9", ls="None", marker=".", zorder=-10)


min_overlap = 0.0001
smap = plt.cm.ScalarMappable(cmap="magma_r")
for x, es, os in zip(x_values, energies, overlaps):
    inds = np.argwhere(os > min_overlap).flatten()
    inds = inds[np.argsort(os[inds])]
    if len(inds) > 0:
        alpha = 1 - np.maximum(np.log(os[inds]) / np.log(min_overlap), 0)
        ax.scatter([x] * len(es[inds]), es[inds], c=os[inds], alpha=alpha, s=15, vmin=0, vmax=1, cmap=smap.cmap)
fig.colorbar(smap, ax=ax, label="Overlap with state of interest")


plt.show()
