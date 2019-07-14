import pytraj as pt
import matplotlib.pyplot as plt
from pyemma.plots import plot_free_energy
import numpy as np
import os

plt.style.use("paper")

# trajfiles = ["sim1/level_01.nc", "sim2/level_01.nc"]
# e2e  = np.empty(0)

# for trajfile in trajfiles:
#     traj = pt.load(trajfile, top="rna_no_water.prmtop")
#     ref = pt.load("rna_h_native.pdb")
#     e2e_i = pt.distance(traj, mask=":3@C1', :6@C1'")
#     e2e = np.hstack((e2e, e2e_i))
#     np.savetxt(os.path.join(os.path.dirname(trajfile), "e2e.txt"), e2e_i)

# #rmsd = pt.rmsd(traj = traj, mask=":3, :4, :5, :6", ref=ref, ref_mask="@C3'")
# ermsdfiles = ["sim1/ermsd.01", "sim2/ermsd.01"]
# ermsd = np.empty(0)
# for efile in ermsdfiles:
#     ermsd = np.hstack((ermsd, np.loadtxt(efile, usecols=1)))
# #radgyr = pt.radgyr(traj=traj, mask="!:Na+")

# print(e2e.shape, ermsd.shape)
# data = np.hstack((ermsd.reshape(-1,1), e2e.reshape(-1,1)))
# np.savetxt("ermsd_e2e.txt", data)
ermsd, e2e = np.loadtxt("./ermsd_e2e.txt", unpack=True)
ax = plot_free_energy(ermsd, e2e, kT=0.6, cbar_label="Free energy (kcal/mol)")
plt.xlabel(r"eRMSD from native")
plt.ylabel(r"C1'/C1' end-to-end distance (\AA)")
plt.tight_layout()
plt.savefig("free_energy.png")
plt.show()
