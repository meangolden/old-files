# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 15:11:15 2021

@author: cp17593
"""
#for the purpose of https://www.overleaf.com/9126177595fdmmkwnfmtvg
#lightweight protocol for RFID, Rayleigh channels
#plotting:  R_{ch}(d) = \frac{\sin(kd - (\theta_{LOS} + \angle{\Gamma} )}{kd}
#Notation:
#a,β: polar and azimuthal angle


import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, BSpline
from scipy import signal
import scipy.special as sp
from scipy.special import jv #ethe Bessel function of the first kind of real order and complex argument
from math import pi
def noise(sigma):
    '''Compute additive white noise of zero mean and variance equal to sigma.
    Returns complex white noise '''

    return np.random.normal(0, sigma/2) + np.random.normal(0, sigma/2)*1j

def channel(sigma):
    '''Compute complex channel coefficient of zero mean and variance equal to sigma. '''
    return np.random.normal(0, sigma/2) + np.random.normal(0, sigma/2)*1j

#fc = 1e9 #carrier frequency in Hz
#l =  3*1e8/fc      #wavelength
l=1 #wavelength

gamma = 1
D = np.linspace(0,6*l,100)
#spatialCorr = [sp.jv(0,2*pi*d/l) for d in D ]
#spatialCorrsqr = [np.sqrt((sp.jv(0,2*pi*d/l))**2) for d in D ]
channelCorrelation_3D = [max(np.sin(2*pi*d/l)/(2*pi*d/l),-np.sin(2*pi*d/l)/(2*pi*d/l)) for d in D]
channelCorrelation_reflection = [max(np.cos(gamma+2*pi*d/l)*np.sin(2*pi*d/l)/(2*pi*d/l),-np.cos(gamma+2*pi*d/l)*np.sin(2*pi*d/l)/(2*pi*d/l))for d in D]
plt.plot(D,channelCorrelation_3D)

plt.plot(D,channelCorrelation_reflection)
#plt.plot(D,mostCommonCorrabs) #absolute value

#plt.plot([l/2,l/2],[-1,1])
#plt.plot([0.38*l,0.38*l],[-1,1])
plt.show()


#CHECK https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1330410
#AND https://www.mendeley.com/reference-manager/reader/b6cae3f3-f127-3c5a-b950-2292fefde4da/0fa75d5f-d80b-67ed-29fd-a8f2e8925847
#AND: Good one:https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1043864
#I would say to validate proximity once the correlation is bigger than 0.35 or so. 
#https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1043864
#https://asa.scitation.org/doi/pdf/10.1121/1.1908122
#https://www.mendeley.com/reference-manager/reader/3ee83584-7d05-3989-bd41-8cbc7298998e/3bcf44e2-b152-6641-361c-a44b22faa433
#^^^^^^ FANTASTIC ^^^^^^^^