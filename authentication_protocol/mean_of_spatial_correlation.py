# -*- coding: utf-8 -*-
"""
Created on Mon May 10 15:49:06 2021

@author: cp17593
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 19:41:51 2021

@author: cp17593
"""
#Everything here is 2D.
#Recreating the graphs of paper "Effect of fading correlation on adaptive arrays in digital mobile radio."
#https://www.mendeley.com/reference-manager/reader/3ee83584-7d05-3989-bd41-8cbc7298998e/3bcf44e2-b152-6641-361c-a44b22faa433
#Empirical evaluation of 2D Rayleigh channels.

import matplotlib.pyplot as plt
import numpy as np
from math import pi
import scipy.special as sp
import statistics as stats

def correlation_sample(D, no_paths, phi, Delta):
    y1=np.zeros(len(D))
    for i in range(len(D)):
        for j in range(no_paths):
            phase = np.random.uniform(phi - Delta, phi + Delta )
            phase_diff = 2*pi*D[i] *np.sin(phase)
            y1[i] += np.cos(phase_diff)
    return y1
        

no_paths, no_samples = 50,7000#
mini_sample = 5
D = np.linspace(0,4,100)   # distance is normalised to wavelenght
Delta = 15 #let's make this narrow
#for the case when theta is in one plane
#phi = pi/2
phi, degree_var = 0, 15
Delta = degree_var * pi/180

#main body
Ysamples=[]
Ysamples =[correlation_sample(D, no_paths, phi, Delta)/no_paths for i in range(no_samples)]
#each row stores a series of samples across the distance (D in size)

Yvar = [stats.variance([Ysamples[r][c] for r in range(no_samples)]) for c in range(len(D))]
Y95percent = [np.percentile([Ysamples[r][c] for r in range(no_samples)],95) for c in range(len(D))] 
Y05percent = [np.percentile([Ysamples[r][c] for r in range(no_samples)],5) for c in range(len(D))] 

l=1
spatialCorr = [sp.jv(0,2*pi*d/l) for d in D ]

if degree_var/90==0:
    plt.plot(D, spatialCorr, label ="Mean value",color="black")
else:
    meanEst = [np.mean([Ysamples[r][c] for r in range(no_samples)]) for c in range(len(D))]
    plt.plot(D, meanEst, label ="Mean value",color="black")
    
SampleMean = [np.mean([Ysamples[r][c] for r in range(mini_sample)]) for c in range(len(D))]
plt.plot(D,Ysamples[0], color="green",label ="One realisation")
plt.plot(D, SampleMean, label ="5-block sample mean",color="red")

plt.plot(D, Y95percent, color="blue", linestyle="--",label =" 5% and 95% percentiles")
plt.plot(D, Y05percent, color="blue", linestyle="--")
plt.legend()
plt.xlabel("$d/\lambda$",fontsize=14)
plt.ylabel("Re($R_{1,2})$",fontsize=14)


