import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

plt.style.use("paper")
params_dict = {
    "xtick.labelsize": 30,
    "ytick.labelsize": 30,
    "axes.labelsize": 30
}

top=0.88,
bottom=0.11,
left=0.555,
right=0.98,
hspace=0.2,
wspace=0.2

if __name__ == "__main__":
    mpl.rcParams.update(params_dict)
    gardient = np.linspace(0, 1, 256).reshape(-1, 1)[::-1]
    gardient = np.hstack((gardient, gardient))

    fig, ax = plt.subplots(figsize=(1, 6))
    ax.imshow(gardient, aspect="auto", cmap=plt.get_cmap("gnuplot_r"))
    ax.xaxis.set_visible(False)
    ax.set_yticks(np.linspace(0, 256, 5))
    ax.set_yticklabels(np.linspace(0, 1, 5)[::-1])
    # plt.tight_layout()
    plt.savefig("colrbar.eps", bbox_inches='tight')
    plt.show()
