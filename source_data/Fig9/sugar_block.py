import numpy as np
import matplotlib.pyplot as plt
import matplotlib

plt.style.use("paper")
ratio = 4

f = plt.figure(figsize=(8,4.5))
gs = plt.GridSpec(ratio, ratio+1) # 4 * 5
ax_ref = f.add_subplot(gs[3, :-1])
ax2_ref = f.add_subplot(gs[3, -1], sharey = ax_ref)
plt.setp(ax2_ref.get_xticklabels(), visible=False)
plt.setp(ax2_ref.get_yticklabels(), visible=False)
plt.setp(ax2_ref.yaxis.get_majorticklines(), visible=False)

if __name__ == "__main__":
    for i in range(3):
        ax = f.add_subplot(gs[i, :-1], sharex = ax_ref)
        ax2 = f.add_subplot(gs[i, -1], sharey = ax, sharex=ax2_ref)
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.setp(ax2.yaxis.get_majorticklines(), visible=False)
        plt.setp(ax2.xaxis.get_majorticklines(), visible=False)

    ax = f.get_axes()
    sugar = np.loadtxt("sugar_angle.txt")
    sugar2 = np.loadtxt("sugar_angle_sim2.txt")
    sugar = np.vstack((sugar, sugar2))
    nsnap = len(sugar[:, 0])
    stride = 50
    x = np.arange(0, nsnap/stride)

    ax[0].scatter(x, sugar[:, 3][::stride], s=4, c='k')
    ax[2].scatter(x, sugar[:, 0][::stride], s=4, c='k')
    ax[4].scatter(x, sugar[:, 1][::stride], s=4, c='k')
    ax[6].scatter(x, sugar[:, 2][::stride], s=4, c='k')

    height, x = np.histogram(sugar[:, 3], 40)
    ax[1].plot(height, (x[1:] + x[:-1]) / 2, color="b")
    ax[1].text(0.7, 0.7, r'\textbf{(d)}', transform=ax[1].transAxes, fontsize=16)

    height, x = np.histogram(sugar[:, 0], 40)
    ax[3].plot(height, (x[1:] + x[:-1]) / 2, color="b")
    ax[3].text(0.7, 0.7, r'\textbf{(a)}', transform=ax[3].transAxes, fontsize=16)

    height, x = np.histogram(sugar[:, 1], 40)
    ax[5].plot(height, (x[1:] + x[:-1]) / 2, color="b")
    ax[5].text(0.7, 0.7, r'\textbf{(b)}', transform=ax[5].transAxes, fontsize=16)

    height, x = np.histogram(sugar[:, 2], 40)
    ax[7].plot(height, (x[1:] + x[:-1]) / 2, color="b")
    ax[7].text(0.7, 0.7, r'\textbf{(c)}', transform=ax[7].transAxes, fontsize=16)
    # ax[1].hist(sugar[:, 3], 40, orientation="horizontal", density=True, color='k')
    # ax[3].hist(sugar[:, 0], 40, orientation="horizontal", density=True, color='k')
    # ax[5].hist(sugar[:, 1], 40, orientation="horizontal", density=True, color='k')
    # ax[7].hist(sugar[:, 2], 40, orientation="horizontal", density=True, color='k')

    for i in range(0, 8, 2):
        ax[i].set_ylim(0, 2*np.pi)
    ax_ref.set_xlim(0, nsnap/stride)
    ax_ref.set_xticks(np.linspace(0, nsnap//stride, 9))
    ax_ref.set_xticklabels(np.arange(9))
    ax_ref.set_xlabel(r"\textbf{Simulation time }($\mu$s)")

    plt.text(0.1, 0.5, 'Pucker Phase(rad)',
            horizontalalignment='right',
            verticalalignment='center',
            rotation='vertical',
            transform=f.transFigure,
            fontsize=20)


    plt.tight_layout()
    f.subplots_adjust(hspace=0)
    f.subplots_adjust(wspace=0)
    plt.savefig("pucker_phase.png", bbox_inches='tight')
    plt.show()
