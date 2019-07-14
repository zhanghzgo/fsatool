import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import networkx as nx
import glob
from networkx.drawing.nx_agraph import to_agraph

def plot_flux_network(pi, fluxmat, fraction=1.0, file="temp.png", direction=True, layout="dot"):
    """ plot the flux network given the pi, tptflux, and fraction"""
    if direction:
        flux = nx.DiGraph()
    else:
        flux = nx.Graph()
    flux.graph["graph"] = {"rankdir": "LR"}
    pi = pi/sum(pi) * len(pi)/2.0


    # set the flux width
    fluxmatwidth = fluxmat/np.sum(fluxmat)*5
    flux = network_add_node(flux, pi)
    pos = nx.spring_layout(flux)
    edge = [(i, j, {
        'penwidth': max(fluxmatwidth[i, j],0.05)*50,
        'label':"{:.2e}".format(fluxmat[i, j]),
        "color":"red", 
        'fontsize':"100", 
        "style":"bold",
        "arrowsize": adjust_arrowsize(fluxmat[i,j])})
        for i in range(len(pi)) for j in range(len(pi)) if fluxmat[i, j] > 0]
    flux.add_edges_from(edge)
    # draw_matplotlib_figure(flux)
    A = to_agraph(flux)
    A.layout(layout)
    A.write("file.dot")
    A.draw(file)

def adjust_arrowsize(fluxmatwidth):
    if(fluxmatwidth > 0.05):
        return 0.05
    else:
        return 4

def network_add_node(flux,pi):
    pngfile = sorted(glob.glob("state*.png"))
    nodelist = []
    paramdict = {}
    for i in range(len(pi)):
        paramdict = {
        'fontsize':"30", 
        'image':pngfile[i],
        "label": "",
        # 'label': "state{:02d}".format(i+1),
        'shape': "none"
        }
        nodelist.append((i, paramdict))
    flux.add_nodes_from(nodelist)
    return flux

if __name__ == "__main__":
    data = np.loadtxt("macro_net_flux.txt")
    nmacro = data.shape[0]
    pi = np.ones(nmacro)/nmacro
    plot_flux_network(pi, data)
