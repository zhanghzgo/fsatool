import numpy as np
import matplotlib.pyplot as plt

plt.style.use("paper")

fig, ax = plt.subplots(1,1, figsize=(8,4), sharex=True)
ermsd01 = np.loadtxt("ermsd.01", usecols=(1))
ermsd02 = np.loadtxt("ermsd.02", usecols=(1))

chi_angle = np.loadtxt("chi_angles")
ermsd = np.hstack((ermsd01, ermsd02))
chi_angle = np.hstack((chi_angle[:,0], chi_angle[:, 1]))

t = np.arange(len(ermsd))
line1, = ax.plot(t[::50], ermsd[::50], label="eRMSD", color="k")
ax.set_ylabel("eRMSD")

ax_sub1 = ax.twinx()
ax_sub1.set_ylabel(r"$\textbf{G}_{\textbf{L4}} \chi$ (rad)", color="r")
line2, = ax_sub1.plot(t[::50], chi_angle[::50], "r--", label=r"$\chi$")
ax_sub1.set_ylim(0, 7)
ax.legend((line1, line2), ("eRMSD", r"$\chi$"), loc='lower right', fontsize="medium", bbox_to_anchor=(0.85, 0))


# ax[1].plot(t[::10], ermsd02[::10], label="Simulation 2", color="k")
# ax[1].set_ylabel("eRMSD")
# ax_sub2 = ax[1].twinx()
# ax_sub2.set_ylabel(r"$\textbf{G}_{\textbf{L4}} \chi$ (rad)", color="r")
# ax_sub2.plot(t[::10], chi_angle[::10, 1], c= "r", label=r"$\chi$")
# ax_sub2.set_ylim(0, 7)

plt.xlim(0, len(ermsd))
plt.xticks(np.linspace(0, len(ermsd), 9), np.arange(9))
ax.set_xlabel(r"Simulation time$(\mu s)$")
plt.tight_layout()
# fig.subplots_adjust(hspace=0)
plt.savefig("ERMSD_chi.png")
plt.show()
