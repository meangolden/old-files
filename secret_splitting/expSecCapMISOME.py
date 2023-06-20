#!/usr/bin/env python3

# Chrys 2020-04-02
# Secrecy capacity against number of antennnas at Eve or Transmit Power. 
# MISOME model. 
# Functions expmisomeSplitCs and expmisomeCs are used from SS_Functions.

# Standard library
from math import *

# PyPI modules
import numpy as np
import seaborn as sns
import matplotlib
#import scipy.linalg as sla
#matplotlib.use("Agg") # Don't require X11.
import matplotlib.pyplot as plt

# Local modules
from SS_Functions import *

#x: Int = foobarbaz()

if __name__ == "__main__":

    #What do you want to plot? A) Sec Capaciry Vs power B) Sec Capcity Vs NE
    method = int(input(" Enter 1 to for a graph against Tx Power\n\
    Enter 2 for a graph against the number of antennas at Eve\n"))
    assert method == 1 or 2

    NA = 4
    xe, ye = .5, 0 #Eve'slocation
    sigma2e1, sigma2e2 = channelStats(xe, ye)
    nIter = 10**2

    if method == 1:

        NE = 3
        PdbInput = np.linspace(-10, 20, 5)

        tradCsVsPower = [expmisomeCs(sigma2e1, 10**(Pdb/10), NA, NE, nIter) for Pdb in PdbInput] # (..) generator expression rather than [...]: list comprehension.
            #N small things better than 1 big thing.

        secsplitCsVsPower = [expmisomeSplitCs(sigma2e1,
                                         sigma2e2,
                                         10**(Pdb/10)/2,
                                         NA//2,
                                         NA//2,
                                         NE, nIter) \
                             for Pdb in PdbInput]


        #noEveCvsPower = [noeveChannCap(10**(Pdb/10), 1) for Pdb in PdbInput]

        fignum = NE
        fig = plt.figure(fignum, figsize=(8, 5)) # figsize dimensions in inches.

        plt.plot(PdbInput, tradCsVsPower, label = 'Conventional')
        plt.plot(PdbInput, secsplitCsVsPower, label = 'Secret Splitting')
        #plt.plot(PdbInput, noEveCvsPower, label = 'no Eve Conventional')
        plt.scatter(PdbInput, tradCsVsPower)
        plt.scatter(PdbInput, secsplitCsVsPower)
        plt.legend()
        plt.title("NE={} , NA={}".format(NE, NA))
        plt.ylabel("Secrecy Capacity")
        plt.xlabel("P (dB)")

        plt.show()
        plt.close()

        #breakpoint()
    elif method == 2:
        P = 10
        NErange = range(1, 16)
        tradCsVsNE = [expmisomeCs(sigma2e1, P, NA, NE, nIter) for NE in NErange]

        secsplitCsVsNE = [expmisomeSplitCs(sigma2e1, sigma2e2, P/2, NA//2, NA//2, NE, nIter) for NE in range(1,16)]

        fignum = 1
        fig = plt.figure(fignum, figsize=(8, 5)) # figsize dimensions in inches.

        plt.plot(NErange, tradCsVsNE, label = 'Conventional')
        plt.plot(NErange, secsplitCsVsNE, label = 'Secret Splitting')
        #plt.plot(PdbInput, noEveCvsPower, label = 'no Eve Conventional')
        plt.scatter(NErange, tradCsVsNE)
        plt.scatter(NErange, secsplitCsVsNE)
        plt.legend()
        plt.title("P={}dB , NA={}".format(10*np.log10(P), NA))
        plt.ylabel("Secrecy Capacity")
        plt.xlabel("NE")

        plt.show()
        plt.close()

    #breakpoint()
