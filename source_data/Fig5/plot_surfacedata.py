import numpy as np
import matplotlib.pyplot as plt
import sys
filename = sys.argv[1]
nbins = int(sys.argv[2])


plt.style.use("paper")
cmap = plt.cm.get_cmap("nipy_spectral")
pot,x,y = np.loadtxt(filename, skiprows=3, usecols=(0,2,3), unpack=True)
x = x.reshape(nbins, nbins)
y = y.reshape(nbins, nbins)
pot = pot.reshape(nbins, nbins)

nonzero = pot.nonzero()
pot[pot==0] = np.inf
pot[nonzero] -= np.min(pot[nonzero])

plt.xlabel(r"Number of Hydrogen Bonds")
plt.ylabel(r"RMSD(\AA)")
plt.ylim(0.1, 20)
plt.contourf(x,y,pot,100,extend="neither", cmap=cmap)
bar = plt.colorbar()
bar.set_label("Free energy (kcal/mol)")
plt.tight_layout()
plt.savefig("free_energy_surface.png")
plt.show()
