#!/usr/bin/env python3

# Chrys 2020-05-11
# something is wrong.
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
    noIter = 10**2
    NA, NAtrad, NB, NE = 2, 4, 1, 4 #both A1 and A2 have NA antennas each.
    phi = 0.8 # phi:(1-phi) is the ratio of Tx power for inf signal : Tx P for AN
    assert NB == 1, "this model works for single antenna Bob"  
    assert NA > 1, "if Alice has only one antenna, she can't generate AN.\
                        you need to ammend your code if you want to include this case."
      
    #Eve's locations
    a = 0.5
    assert a != 1 #do not place Eve right on the Base Stations. 
    evesLoc = [[a, 0], [-a, 0]]  #the eavesdroppers' locations
    K = len(evesLoc) # the number of eavesdroppers. 
    
    Pdb = np.linspace(-20, 100, 30)
    P = 10**(Pdb/10)
    nP = len(P)
    #initialisation
    CsSplit = [np.zeros(nP)] * noIter
    CsTrad = [np.zeros(nP)] * noIter
    CsTradbound = [np.zeros(nP)] * noIter
    
    #Channel statistics of the eavesdroppers.
    channelsStats = np.empty([K, 2]) 
    for k in range(K):
        channelsStats[k][0], channelsStats[k][1] = SS.channelStats(evesLoc[k][0], evesLoc[k][1])
    EvesforA1 = []
    EvesforA2 = []
    #Divide the eavesdroppers in two subsets according to their channel statistics. 
    for k in range(K):
        if channelsStats[k][0] < channelsStats[k][1]:
            EvesforA1.append(channelsStats[k][0])
        else:
            EvesforA2.append(channelsStats[k][1])
    #Find the worst channel statistic for A1, repeat for A2.        
    maxS2EA1, maxS2EA2  = max(EvesforA1), max(EvesforA2)
        
    for iter in range(noIter):
        #Bob's channels
        Hb1, Hb2 = SS.channelMatrix(NB, NA, 1), SS.channelMatrix(NB, NA, 1)
        Hb1trad = SS.channelMatrix(NB, NAtrad, 1)
        #Eve's channels
        He1, He2 = SS.channelMatrix(NE, NA, maxS2EA1), SS.channelMatrix(NE, NA, maxS2EA2)
        #breakpoint()
        He1trad = SS.channelMatrix(NE, NAtrad, max(channelsStats.T[0]))
        #I keep the first Base Station (no optimal single BS choice.)
        #print("worst for A1", maxS2EA1)
        #print("worst for A2", maxS2EA2)
        #print("worst for Atrad", max(channelsStats.T[0]))
        
        #Beamformers and noisevectors
        maxEigVal1, beamformer1, noiseVctr1, ANcov1 = SS.transmittingwithAN(Hb1)
        maxEigVal2, beamformer2, noiseVctr2, ANcov2 = SS.transmittingwithAN(Hb2)
        maxEigValtrad, beamformertrad, noiseVtrad, ANcovtrad = SS.transmittingwithAN(Hb1trad)
        
        #SNRb1 = [SS.SNRBob1(phi, p/2, Hb1) for p in P]
        #SNRb2 = [SS.SNRBob1(phi, p/2, Hb1) for p in P]
        #SNRbtrad = [SS.SNRBob1(phi, p, Hb1trad) for p in P]
        #or
        
        SNRb1 = [SS.SNRBob(phi, p/2, maxEigVal1) for p in P]
        SNRb2 = [SS.SNRBob(phi, p/2, maxEigVal2) for p in P]
        SNRbtrad = [SS.SNRBob1(phi, p, maxEigValtrad) for p in P]
        #print("SNRsBob", SNRb1, SNRb2, SNRbtrad)
       
        #SNRe1 = [SS.SNREve1(phi, p/2, beamformer1, noiseVctr1, He1) for p in P]
        #SNRe2 = [SS.SNREve1(phi, p/2, beamformer2, noiseVctr2, He2) for p in P]
        #SNRetrad = [SS.SNREve1(phi, p, beamformertrad, noiseVtrad, He1trad) for p in P]
        #or
        SNRe1 = [SS.SNREve(phi, p/2, beamformer1, ANcov1, He1) for p in P]
        SNRe2 = [SS.SNREve(phi, p/2, beamformer2, ANcov2, He2) for p in P]
        SNRetrad = [SS.SNREve(phi, p, beamformertrad, ANcovtrad, He1trad) for p in P]
        #print("SNRsEve", SNRe1, SNRe2, SNRetrad)
        #capacities
        CB1, CB2 = [log2(1 + b1) for b1 in SNRb1], [log2(1 + b2) for b2 in SNRb2]
        CE1, CE2 = [log2(1 + e1) for e1 in SNRe1], [log2(1 + e2) for e2 in SNRe2]
        C1 = [max(0, cb1 - ce1) for (cb1, ce1) in zip(CB1, CE1)]#, max(0, CB2 - CE2)
        C2 = [max(0, cb2 - ce2) for (cb2, ce2) in zip(CB2, CE2)]
        
        #Secrecy rates
        CsSplit[iter] = [SS.splitCs(c1, c2) for (c1, c2) in zip(C1, C2)]
        CsTrad[iter] = [max(0, log2(1 + b) - log2(1 + e)) for (b, e) in zip(SNRbtrad, SNRetrad)]
        #print("number of iteration:", iter)
        ##print("CsSplit=\kn", np.round(CsSplit,3))
        #print("CsTrad=\n", np.round(CsTrad,3))
        
        #breakpoint()
        
        #boundary for Ctrad (although secrecy rate is more appropriate.)
        CsTradbound[iter] = [SS.misomeCs(Hb1trad, He1trad, p) for p in P]
        
    #Expected Secrecy Rates
    ##breakpoint()
    expCsSplit = [np.mean(a) for a in zip(*CsSplit)]
    expCsTrad = [np.mean(a) for a in zip(*CsTrad)]
    expCsboundTrad = [np.mean(a) for a in zip(*CsTradbound)] 
    

    #figure
    fignum = NE
    fig = plt.figure(fignum, figsize=(8, 5)) # figsize dimensions in inches.

    plt.plot(Pdb, expCsSplit, label = 'Secret Splitting')
    plt.plot(Pdb, expCsTrad, label = 'Conventional secrecy coding')
    plt.plot(Pdb, expCsboundTrad, label = 'Bound for Conven.')
    plt.legend()
    plt.title("NE={}, eav. loc={}, NA1=NA2={}, NAtrad={}, phi={}".format(NE, evesLoc, NA, NAtrad, phi))
    plt.ylabel("Expected Secrecy Rate")
    plt.xlabel("total tx power (dB)")

    plt.show()
    plt.close() 
    plt.savefig("NE{}NA{}phi{}.pdf".format(NE,NA, np.round(phi*10,1)))
    ##breakpoint() #you may want to #print the channel realisations of the eavesdroppers. 



