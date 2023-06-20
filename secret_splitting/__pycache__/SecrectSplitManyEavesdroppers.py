#!/usr/bin/env python3

# Chrys 2020-05-04
# MISOME channel. misomeCs and splitsCs functions from SS_Functions are used.
# Deterministic position of the eavesdroppers. >> Line 40
# Bob is positioned at O(0,0). Base stations at A1(1,0) and A2(-1,0).

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
#this code performs well!
    noIter = 10**2
    NA, NAtrad, NB, NE = 2, 4, 1, 4 #both A1 and A2 have NA antennas each.
    assert NB == 1, "this model works for single antenna Bob"
    
    Pdb = np.linspace(-10, 100, 20)
    P = 10**(Pdb/10)
    nP = len(P)
    CsSplit_ = [0] * nP
    CsTrad_ = [0] * nP
    
    #location and channel stats of the eavesdroppers.
    a = 0.3
    evesLoc = [[a, 0], [-a, 0]] #the eavesdroppers' locations
    K = len(evesLoc) # the number of eavesdroppers. 
    channelsStats = np.empty([K, 2]) # the kth element of this array stores the
    
    for k in range(K):
        channelsStats[k][0], channelsStats[k][1] = SS.channelStats(evesLoc[k][0], evesLoc[k][1])    
        
    for iter in range(noIter):
   
        Hb1, Hb2 = SS.channelMatrix(1, NA, 1), SS.channelMatrix(1, NA, 1)
        Hb1trad = SS.channelMatrix(1, NAtrad, 1)
        
        He1s = [SS.channelMatrix(NE, NA, channelsStats[k][0]) for k in range(K)]
        He2s = [SS.channelMatrix(NE, NA, channelsStats[k][1]) for k in range(K)]
        LA1s, LA2s = [LA.norm(He1s[k]) for k in range(K)], [LA.norm(He2s[k]) for k in range(K)]
        
        He1trad = [SS.channelMatrix(NE, NAtrad, channelsStats[k][0]) for k in range(K)]
        LA1trad = [LA.norm(He1trad[k]) for k in range(K)]
        worstchannelA1trad = He1trad[np.argmax(LA1trad)]
        
        I1 = [int(LA1s[k] < LA2s[k]) for k in range(K)]
        I2 = (I1 + np.ones(K))%2

        LA1s, LA2s = [LA1s[k]*I1[k] for k in range(K)], [LA2s[k]*I2[k] for k in range(K)]
        worstChannelforA1 = He1s[np.argmax(LA1s)]
        worstChannelforA2 = He2s[np.argmax(LA2s)]
    
        if np.sum(I1) == 0:
            worstChannelforA1 = [] 
        if np.sum(I2) == 0:
            worstChannelforA2 = []
        #up to this point we attained the appropriate channel matrices

        C1 = [SS.misomeCs(Hb1, worstChannelforA1, p/2) for p in P]
        C2 = [SS.misomeCs(Hb2, worstChannelforA2, p/2) for p in P]
        CsSplit = [SS.splitCs(C1[i], C2[i]) for i in range(nP)]
        CsTrad = [SS.misomeCs(Hb1trad, worstchannelA1trad, p) for p in P]
        CsSplit_ = [CsSplit_[i] + CsSplit[i] for i in range(nP)]
        CsTrad_ = [CsTrad_[i] + CsTrad[i] for i in range(nP)]

        
    CsSplits = [i/noIter for i in CsSplit_]
    CsTrads = [i/noIter for i in CsTrad_]

    #figure
    fignum = NE
    fig = plt.figure(fignum, figsize=(8, 5)) # figsize dimensions in inches.

    plt.plot(Pdb, CsSplits, label = 'Secret Splitting')
    plt.plot(Pdb, CsTrads, label = 'Conventional secrecy coding')
    plt.legend()
    plt.title("NE={}, eav. loc={}, NA1=NA2={}, NA={}".format(NE, evesLoc, NA, NAtrad))
    plt.ylabel("Expected Secrecy Rate")
    plt.xlabel("total tx power (dB)")
    plt.show()
    plt.close()
    
    #breakpoint()
    
   
    
 


