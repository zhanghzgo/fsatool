import numpy as np
import matplotlib.pyplot as plt
plt.style.use('paper')

plt.subplots(1,1,figsize=(6,4))

data = np.loadtxt("check.out", skiprows=1, comments="#")[:, :]
font = {"family": "Times New Roman"}
lw = 2; mz = 4

if __name__ == "__main__":
    plt.plot(data[:, 0]*10, data[:, 1]*10, 'ko-', linewidth=lw, markersize=mz)
    plt.plot(data[:, 0]*10, data[:, 2]*10, 'ko-', linewidth=lw, markersize=mz)
    plt.plot(data[:, 0]*10, data[:, 3]*10, 'ko-', linewidth=lw, markersize=mz)
    plt.plot(data[:, 0]*10, data[:, 4]*10, 'ko-', linewidth=lw, markersize=mz)
    plt.yscale('log')
    plt.ylim(20,100)
    plt.xlim(0, 400)
    plt.xlabel("lag time (ps)", fontdict=font)
    plt.ylabel("implied time scale(ps)", fontdict=font)
    plt.tight_layout()
    plt.savefig("timescale.png")
    plt.show()
