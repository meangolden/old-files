# -*- coding: utf-8 -*-
"""
Created on Tue May 18 21:01:23 2021

@author: cp17593
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 17:52:57 2021

@author: cp17593
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, BSpline
from scipy import signal
import scipy.special as sp
from scipy.special import jv #ethe Bessel function of the first kind of real order and complex argument
from math import pi



def correlation_sample(D, no_paths, phi, Delta):
    y1=np.zeros(len(D))
    for i in range(len(D)):
        for j in range(no_paths):
            phase = np.random.uniform(phi - Delta, phi + Delta )
            phase_diff = 2*pi*D[i] *np.sin(phase)
            y1[i] += np.cos(phase_diff)
    return y1
        

no_paths, no_samples = 50,8000#
D = np.linspace(0,4,100)   # distance is normalised to wavelenght

#for the case when theta is in one plane
#phi = pi/2
phi, degree_var = 0, 30
Delta = degree_var * pi/180

#main body
Ysamples=[]
Ysamples =[correlation_sample(D, no_paths, phi, Delta)/no_paths for i in range(no_samples)]
#each row stores a series of samples across the distance (D in size)

meanEst = [np.mean([Ysamples[r][c] for r in range(no_samples)]) for c in range(len(D))]

mostCommonCorr = [np.sin(2*pi*d)/(2*pi*d) for d in D]

plt.plot(D, mostCommonCorr, label ="without dominant path")

plt.xlabel("$d/\lambda$",fontsize=12)
plt.ylabel("Channel correlation",fontsize=12)
plt.plot(D,meanEst,label ="with dominant path")
plt.legend(fontsize=12)

plt.show()



#TO DO NEXT:
#CHECK https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1330410
#AND https://www.mendeley.com/reference-manager/reader/b6cae3f3-f127-3c5a-b950-2292fefde4da/0fa75d5f-d80b-67ed-29fd-a8f2e8925847
#AND: Good one:https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1043864
#I would say to validate proximity once the correlation is bigger than 0.35 or so. 
#https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1043864
#https://asa.scitation.org/doi/pdf/10.1121/1.1908122
#https://www.mendeley.com/reference-manager/reader/3ee83584-7d05-3989-bd41-8cbc7298998e/3bcf44e2-b152-6641-361c-a44b22faa433
#^^^^^^ FANTASTIC ^^^^^^^^