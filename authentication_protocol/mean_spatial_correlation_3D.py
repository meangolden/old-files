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
#Recreating the graphs of paper "Effect of fading correlation on adaptive arrays in digital mobile radio."
#https://www.mendeley.com/reference-manager/reader/3ee83584-7d05-3989-bd41-8cbc7298998e/3bcf44e2-b152-6641-361c-a44b22faa433
#Empirical evaluation of 3D Rayleigh channels.

import matplotlib.pyplot as plt
import numpy as np
from math import pi
import scipy.special as sp
import statistics as stats

def correlation_sample(D, no_paths):
    y2=np.zeros(len(D))
    for i in range(len(D)):
        for j in range(no_paths):
            phi= np.random.uniform(0,pi) #angle form north
            #theta = np.random.uniform(0,2*pi) #azimuthal angle
            #theta = np.random.uniform(0,0) #azimuthal angle
            phase_difference= 2*pi*D[i]*np.cos(phi)
            y2[i] += (np.pi/2)*np.cos(phase_difference)*np.sin(phi)
    return y2
        



no_paths, no_samples = 50,1000#
mini_sample = 4
D = np.linspace(0,10,100)   # distance is normalised to wavelenght

#for the case when theta is in one plane
#phi = pi/2
phi, degree_var = 0, 90
Delta = degree_var * pi/180

#main body
Ysamples=[]
Ysamples =[correlation_sample(D, no_paths)/no_paths for i in range(no_samples)]
#each row stores a series of samples across the distance (D in size)

Yvar = [stats.variance([Ysamples[r][c] for r in range(no_samples)]) for c in range(len(D))]
Y95percent = [np.percentile([Ysamples[r][c] for r in range(no_samples)],95) for c in range(len(D))] 
Y05percent = [np.percentile([Ysamples[r][c] for r in range(no_samples)],5) for c in range(len(D))] 

l=1
spatialCorr = [np.sin(2*pi*d)/(2*pi*d) for d in D]


plt.plot(D, spatialCorr, label ="$R$",color="black")
SampleMean = [np.mean([Ysamples[r][c] for r in range(mini_sample)]) for c in range(len(D))]
plt.plot(D,Ysamples[0], color="green",label ="$M=1$")
plt.plot(D, SampleMean, label ="$M=4$",color="red")

plt.plot(D, Y95percent, color="blue", linestyle="--",label =" 5% and 95% percentiles")
plt.plot(D, Y05percent, color="blue", linestyle="--")
plt.plot(D, spatialCorr,color="black")
plt.legend()
plt.xlabel("$d/\lambda$",fontsize=12)
plt.ylabel("$\hat{R}$",fontsize=12)


