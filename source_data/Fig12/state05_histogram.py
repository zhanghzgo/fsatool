import barnaba as bb
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

plt.style.use("paper")
params_dict = {
    "xtick.labelsize": 30,
    "ytick.labelsize": 30,
    "axes.labelsize": 30
}
matplotlib.rcParams.update(params_dict)

if __name__ == "__main__":
    angles,res = bb.backbone_angles("state05/state_05.nc", "rna_no_water.prmtop", residues=['G_6_0'], angles=['chi'])
    angles[np.where(angles<0.0)] += 2*np.pi
    plt.xlabel(r"$\textbf{G}_{\textbf{L4}} \quad \chi$ \textbf{(rad)}")
    plt.ylabel(r"$\textbf{Probability density}$")
    plt.ylim(0,1.4)
    plt.annotate('syn', xy=(1.5, 0.6), xytext=(2, 0.8),
                arrowprops=dict(facecolor='black', shrink=0.05),
                 size=35
                )
    # plt.annotate('high-anti', xy=(4.9, 0.43), xytext=(3.7, 0.6),
    #             arrowprops=dict(facecolor='black', shrink=0.05),
    #              size="xx-large"
    #             )

    h, x = np.histogram(angles.ravel(), bins=40, density=True)
    plt.plot((x[:-1] + x[1:])/2, h, 'k')
    plt.tight_layout()
    plt.savefig("state05.eps")
    plt.show()
