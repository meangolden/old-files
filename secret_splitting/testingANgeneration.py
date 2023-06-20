#!/usr/bin/env python3

# Chrys 2020-05-20

# Standard library
from math import *

# PyPI modules
import numpy as np
import seaborn as sns
import matplotlib
#import scipy.linalg as sla
#matplotlib.use("Agg") # Don't require X11.
import matplotlib.pyplot as plt
from numpy import linalg as LA

# Local modules
import SS_Functions as SS

#x: Int = foobarbaz()

if __name__ == "__main__":
    #we know the channel statistics of the eavesdropper only. There is something
    #wrong whith this code. To see why, pick NE = 4 and NAtrad = 4.
    
    #parameters

    NA, NE = 4, 4
    phi = 0.6 # phi:(1-phi) is the ratio of Tx power for inf signal : Tx P for AN
    sigmaE = 1
    #channels    
    Hb = np.random.normal(0,1,(1,NA))
    He = np.random.normal(0,sigmaE,(NE,NA))
    
    #power
    Pdb = np.linspace(-20, 100, 20)
    P = 10**(Pdb/10)
    nP = len(P)
    
    CsTradperfectCSI = [SS.misomeCs(Hb, He, p) for p in P]
    

    CsTradwithAN = [SS.RsANbmfEve(phi, p, Hb, He) for p in P]
    #CsTradwithANNe1 = [SS.RsANbmNe1(phi, p, Hb, He) for p in P]

    
    
    #figure
    fignum = NE
    fig = plt.figure(fignum, figsize=(8, 5)) # figsize dimensions in inches.
    #ax = plt.axes(projection='3d')

    plt.plot(Pdb, CsTradperfectCSI, label = 'perfect CSI conv')
    plt.plot(Pdb, CsTradwithAN, label = 'conv with AN')
    #plt.plot(Pdb, CsTradwithANNe1, label = 'conv with AN, Eve no rx beamforming')
    plt.legend()
    #plt.title("NE={}, eav. loc={}, NA1=NA2={}, NAtrad={}, phi={}".format(NE, evesLoc, NA, NAtrad, phi))
    plt.ylabel("Secrecy Rate")
    plt.xlabel("total tx power (dB)")

    plt.show()
    plt.close() 
    #breakpoint()



