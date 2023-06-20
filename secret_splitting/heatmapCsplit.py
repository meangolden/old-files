#!/usr/bin/env python3

# Chrys 2020-04-02

# Heatmap of the capacity under SS in a MISOME model. Gaussian channels.
# Applies for the case when we do not know Eve's channel.
#Redo everything. How do you define capacity?

from math import *
import numpy as np
import seaborn as sns
import matplotlib
from itertools import product
from SS_Functions import *

#matplotlib.use("Agg") # Don't require X11.
import matplotlib.pyplot as plt

#x: Int = foobarbaz()

if __name__ == "__main__":

    NA1, NA2, Pdb = 4, 4, 10
    NE = 2
    P = 10**(Pdb/10)
    #Place Eve at a meshgrid
    resol = 100
    Xe, Ye = np.linspace(-1.5, 1.5, resol),  np.linspace(-1.5, 1.5, resol)
    SecCapSample = np.empty((len(Xe), len(Ye)))

            
    A = [channelStats(xe, ye) for ye in Ye for xe in Xe]
    SNRB = 1  
    alpha = 2 # pathloss exponent
    Cap1 = [max(0, log2(1 + SNRB) - log2(1 + a)) for a,b in A]
    Cap2 = [max(0, log2(1 + SNRB) - log2(1 + b)) for a,b in A]
    CapSS = [splitCs(C1,C2) for C1,C2 in zip(Cap1,Cap2)]
    C = np.reshape(CapSS, (resol, resol))
    
    # figure
    fig, ax = plt.subplots()
    # resetting the axes to the real values
    dx = (Xe[1]-Xe[0])/2.
    dy = (Ye[1]-Ye[0])/2.
    extent = [Xe[0]-dx, Xe[-1]+dx, Ye[0]-dy, Ye[-1]+dy]
    im = plt.imshow(C, extent=extent)
    cbar = ax.figure.colorbar(im)
    plt.xlabel("x"); plt.ylabel("y")
    plt.text(1, 0,"A1", fontsize=12)
    plt.scatter([0, -1, 1], [0, 0, 0], s = 5000, c=(0,0,0), alpha=0.1)
    plt.text(-1, 0,"A2", fontsize=12)
    plt.text(0, 0,"B", fontsize=12)
    plt.show()




