# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 11:25:13 2021

@author: cp17593
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, BSpline
from scipy import signal

def cosine_similarity(X, Y):
    '''Compute cosine similarity between X and Y.
    Cosine similarity, or the cosine kernel, computes similarity as the
    normalized dot product of X and Y:
        K(X, Y) = <X, Y> / (||X||*||Y||)]
        X : complex number
        Y : complex number'''
    num = X.real*Y.real + X.imag*Y.imag
    den = np.sqrt(X.real**2 + X.imag**2) * np.sqrt(Y.real**2 + Y.imag**2)
    return num/den

def eucl_dist(X, Y):
    '''Compute the Euclidean distance between  X and Y.
        X : complex number
        Y : complex number
        Output: real positive number'''
    s = (X.real - Y.real)**2 + (X.imag - Y.imag)**2
    return np.sqrt(s)

def rmsd(X, Y):
    ''' Compute the root mean square deviation of the vectors of X, Y.
        Returns a positive real number'''
    assert len(X) == len(Y)
    n = len(X)
    diff = [(X[i] - Y[i])**2 for i in range(n)]
    return np.sqrt(sum(diff)/n)
        
    

def noise(sigma):
    '''Compute additive white noise of zero mean and variance equal to sigma.
    Returns complex white noise '''

    return np.random.normal(0, sigma/2) + np.random.normal(0, sigma/2)*1j

def channel(sigma):
    '''Compute complex channel coefficient of zero mean and variance equal to sigma. '''
    return np.random.normal(0, sigma/2) + np.random.normal(0, sigma/2)*1j

n = 2# the number of channel realisations where the transmission took place.
S_leg = []  #cosine similariry
D_leg = []  #euclidena distance
S_attack = []
D_attack = []
corr_colocated = 0  #correlation between a colocated device
corr_attack = 0  #correlation between a distanced device
samples_Bob = []
samples_Alice = []
samples_Eve = []
B, A, E = [], [], []

for i in range(n):
    h, g = channel(1), channel(1)
    b = h + noise(0.5)   #Bob's channel
    B.append(np.sqrt((b.real)**2+(b.imag)**2)) #for RSS measurements
    a = h + noise(0.2)   #when Alice is colocated
    A.append(np.sqrt((a.real)**2+(a.imag)**2)) #for RSS measurements
    e = g + noise(0.2)   #distance
    E.append(np.sqrt((e.real)**2+(e.imag)**2)) #for RSS measurements
    S_leg.append(cosine_similarity(b,a))
    D_leg.append(eucl_dist(b,a))
    S_attack.append(cosine_similarity(b,e))
    D_attack.append(eucl_dist(b,e))
    samples_Bob.append(b)
    samples_Alice.append(a)
    samples_Eve.append(e)
    
corr_colocated = signal.correlate(samples_Bob, samples_Alice, mode='same')
corr_attack= signal.correlate(samples_Bob, samples_Eve, mode='same')
auto_corr = signal.correlate(samples_Bob, samples_Bob, mode='same')  
#print("corr colated =", corr_colocated)
#print("coor attack=", corr_attack)
#print("auto correlation", auto_corr)
#plt.plot(corr_attack)
#plt.plot(corr_colocated)

corr_coeff_auto = np.real(np.corrcoef(samples_Bob, samples_Bob))[0][0]
corr_coeff_colocated = np.real(np.corrcoef(samples_Bob, samples_Alice))[0][1]
corr_coeff_attack = np.real(np.corrcoef(samples_Bob, samples_Eve))[0][1]


#print(rmsd(S_leg, np.ones(n)), rmsd(S_attack, np.ones(n)))
#print(rmsd(D_leg, np.zeros(n)), rmsd(D_attack, np.zeros(n)))
print(round(corr_coeff_auto,3),round(corr_coeff_colocated,3),round(corr_coeff_attack,3))    
#plt.subplot(1, 2, 1)
#plt.scatter(range(n),S_leg)
#plt.subplot(1, 2, 2)
#plt.scatter(range(n), S_attack)
#plt.savefig('cosine_similarity.svg') 
#plt.scatter(range(n),D_leg)
#plt.scatter(range(n), D_attack)
#plt.scatter(range(n),S_leg)
#plt.scatter(range(n), S_attack)
#plt.show()

##RSS measurements
#xnew = np.linspace(0, n-1, 300) 
#spl1 = make_interp_spline(range(n), B, k=3)  # type: BSpline
#plt.axhline(y=0.5, color='black', linestyle='-')
#plt.axhline(y=0.75, color='black', linestyle='-')
#power_smooth1 = spl1(xnew)
#spl2 = make_interp_spline(range(n), A, k=3)  # type: BSpline
#power_smooth2 = spl2(xnew)
#spl3 = make_interp_spline(range(n), E, k=3)  # type: BSpline
#power_smooth3 = spl3(xnew)
#plt.plot(xnew, power_smooth1,color='black')
#ax = plt.gca()
#ax.axes.xaxis.set_ticks([])
#ax.axes.yaxis.set_ticks([])
#plt.savefig('RSSBob.svg') 
#plt.show()
#plt.plot(xnew, power_smooth2)
#plt.axhline(y=0.5, color='black', linestyle='-')
#plt.axhline(y=0.75, color='black', linestyle='-')
#ax = plt.gca()
#ax.axes.xaxis.set_ticks([])
#ax.axes.yaxis.set_ticks([])
#plt.savefig('RSSAlice.svg') 
#plt.show()
#plt.plot(xnew, power_smooth3,color='red')
#plt.axhline(y=0.5, color='black', linestyle='-')
#plt.axhline(y=0.75, color='black', linestyle='-')
#ax = plt.gca()
#ax.axes.xaxis.set_ticks([])
#ax.axes.yaxis.set_ticks([])
#plt.savefig('RSSEve.svg') 
#plt.show()


###############

