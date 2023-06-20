#!/usr/bin/env python3

# Chrys 2020-04-02
#heatmap of the difference

from math import *
import numpy as np
import seaborn as sns
import matplotlib
import itertools
#import scipy.linalg as sla
#matplotlib.use("Agg") # Don't require X11.
import matplotlib.pyplot as plt
from SS_Functions import *

#x: Int = foobarbaz()

def diffCs(CsTrad, CsSplit):
    ''' Returns the difference and the sign of the diffence between the
        conventional Secrecy Capacity and
        Secrecy Capacity under Secret Splitting. Negative sign means Conventional
        scheme perform better. A positive sign indicates that secret splitting performs better.'''
    difference = CsSplit - CsTrad
    assert np.isscalar(difference), (type(difference), difference)
    if difference > 0:
        sign = 1
    else:
        sign = 0
    return difference, sign

if __name__ == "__main__":

    xe, ye = 0.5, 0 #Eve'slocation
    NA_max, NE_max = 16, 16
    nIter = 1
    P = 10 * log10(200) #linear
    sigma2e1, sigma2e2 = channelStats(xe, ye)
   
    #Xe, Ye = np.linspace(-1.5, 1.5, 15), -np.linspace(-1.5, 1.5, 15)
    #Dave's way:
    #minNA, minNE = 2, 1
    #maxNA, maxNE = 16, 16
    #rngNA, rngNE = \
        #list(range(minNA, maxNA+1, 2)), \
       # list(range(minNE, maxNE+1))

    #diff = np.empty((len(rngNA), len(rngNE)))
    #sign = np.empty((len(rngNA), len(rngNE)))

    #for (i,na),(j,ne) in itertools.product(enumerate(rngNA),
    #                                       enumerate(rngNE)):
        #Cs_Trad = expmisomeCs(sigma2e1, P, na, ne, nIter)
        #Cs_Split = expmisomeSplitCs(sigma2e1, sigma2e2, 0, na//2, na//2, ne, nIter)
        
        #diff[i][j], sign[i][j] = diffCs(Cs_Trad, Cs_Split)
    
    #My way
    diff = np.zeros((NA_max // 2, NE_max))
    sign = np.zeros((NA_max // 2, NE_max))
    Cs_Trad_Ma = np.zeros((NA_max//2, NE_max))
    Cs_SS_Ma = np.zeros((NA_max//2, NE_max))
    NA = [2 * (x + 1) for x in range(NA_max // 2)]
    NE =  [x + 1 for x in range(NE_max)]
    
    for i in range(NA_max//2):
        for j in range(NE_max):
            Cs_Trad = expmisomeCs(sigma2e1, P, NA[i], NE[j], nIter)
            Cs_Split = expmisomeSplitCs(sigma2e1, sigma2e2, P / 2, NA[i] // 2, NA[i] // 2, NE[j], nIter)
            diff[i][j], sign[i][j] = diffCs(Cs_Trad, Cs_Split)
            
            #Cs_Trad_Ma[i][j] = Cs_Trad #used for validation of the code only
            #Cs_SS_Ma[i][j] = Cs_Split #same as above
         


    fig, ax = plt.subplots()
    #im = ax.imshow(diff)
    #im = ax.imshow(sign)
    im = ax.imshow(diff)
    ax.set_ylim(0, NA_max//2 - 1)
    #ax.set_xlim(0, NE_max)
    ax.set_xlabel("NE")
    ax.set_ylabel("NA")

    ax.set_xticks(range(NE_max))
    ax.set_yticks(range(NA_max//2))
    ax.set_xticklabels(NE)
    ax.set_yticklabels(NA)
    ax.set_title("Tx Power = %d (linear)" % P)
    plt.show()
    #breakpoint()


