#!/usr/bin/env python3

# Chrys 2020-04-21

# Standard library
from math import *

# PyPI (Python Package Index) library
import numpy as np 
from numpy import linalg as LA
import matplotlib
#matplotlib.use("Agg") # Don't require X11.
import matplotlib.pyplot as plt

# Local library
import SS_Functions as SS

#x: Int = foobarbaz()
#all algorithms have been verified!

if __name__ == "__main__":
    
    method = input(" Press 1 to run the code for one matrix\n Press 2 to get the expected secrecy capacity\n")
    assert int(method) == 1 or 2
    
    if int(method) == 1:

        Hb = SS.channelMatrix(1, 2, 1)
        He = SS.channelMatrix(1, 2, 1)
           
        Pdb = np.linspace(-10, 100, 100)
        
        CsMISOME = [SS.misomeCs(Hb, He, 10**(Pdb/10)) for Pdb in Pdb]
        CsNA2 = [SS.NaTwoCs(Hb, He, 10**(Pdb/10)) for Pdb in Pdb]
        CsMISOSE = [SS.misoseCs(Hb, He, 10**(Pdb/10)) for Pdb in Pdb]
        
        fignum = 1
        fig = plt.figure(fignum, figsize=(8, 5)) # figsize dimensions in inches.

        plt.plot(Pdb, CsMISOME , label = 'method MISOME')
        plt.plot(Pdb, CsNA2, label = '2 tx antennas')
        plt.plot(Pdb, CsMISOSE, label = 'method MISOSE')
        plt.legend(['method MISOME', 'method 2 tx antennas' ,'method MISOSE'])
        plt.title("Verifying algorithms for NA = 2, NB = NE = 1")
        plt.ylabel("Secrecy Capacity")
        plt.xlabel("P (dB)")
        plt.show()
    
    else:
        #Expected Secrecy Capacity
        nIter = 10**3
        NA, NE = 2, 1
        Pdb = np.linspace(-10, 100, 10)
        Cs1_, Cs2_, Cs3_ = [0]*len(Pdb), [0]*len(Pdb), [0]*len(Pdb) #initialisation
        for i in range(len(Pdb)):
            for j in range(nIter):
                Hb = SS.channelMatrix(1, NA, 1) #Bob is a single antenna node
                He = SS.channelMatrix(NE, NA, 1)
                c1 = SS.misomeCs(Hb, He, 10**(Pdb[i]/10))
                c2 = SS.NaTwoCs(Hb, He, 10**(Pdb[i]/10))
                c3 = SS.misoseCs(Hb, He, 10**(Pdb[i]/10))
                #cb = np.log2(1 + P*np.linalg.norm(Hb)**2) #when there is no eavesdropper MRT maximises Bob's SNR
                Cs1_[i] = Cs1_[i] + max(0 , c1)
                Cs2_[i] = Cs2_[i] + max(0 , c2)
                Cs3_[i] = Cs3_[i] + max(0 , c3)
            Cs1_[i], Cs2_[i], Cs3_[i] = Cs1_[i]/nIter, Cs2_[i]/nIter, Cs3_[i]/nIter
        
        breakpoint()
        fignum = 2
        fig = plt.figure(fignum, figsize=(8, 5)) # figsize dimensions in inches.
       
        plt.plot(Pdb, Cs1_, label = 'method MISOME')
        plt.plot(Pdb, Cs2_, label = '2 tx antennas')
        plt.plot(Pdb, Cs3_, label = 'method MISOSE')
        plt.legend(['method MISOME', 'method 2 tx antennas' ,'method MISOSE'])
        plt.title("Verifying algorithms for NA = 2, NB = NE = 1")
        plt.ylabel("Expected Secrecy Capacity")
        plt.xlabel("P (dB)")
        plt.show()
        breakpoint()
