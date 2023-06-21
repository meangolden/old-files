# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 17:52:57 2021

@author: cp17593
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
#from scipy.interpolate import make_interp_spline, BSpline
#from scipy import signal
import scipy.special as sp
#from scipy.special import jv #ethe Bessel function of the first kind of real order and complex argument
from math import pi


###############################################################################
#
figSizeXcm = 12
figSizeYcm = 7

# Font size in Figures
fontSize = 11
#    User Inputs
#
###############################################################################

# Figure Width (cm)

# Line width in Figures (pt)
lineWidth = 1.0

# Choice of Font ("Arial" / "Times New Roman" / "CMU Serif")
fontChoice = "CMU Serif"

# Output File Paths for Plot
outputFilePathPDF = 'figureSizeExample.pdf'

l=1
D = np.linspace(0.001,14*l,10000)

#fc = 1e9 #carrier frequency in Hz
#l =  3*1e8/fc      #wavelength


###############################################################################
#
#    Main Code
#
###############################################################################
spatialCorr = [abs(sp.jv(0,2*pi*d/l)) for d in D ]
mostCommonCorr = [abs(np.sin(2*pi*d)/(2*pi*d)) for d in D]

   
###############################################################################
#
#    Plotting parameters
#
###############################################################################
# Rebuild the matplotlib font cache

matplotlib.font_manager._rebuild()
   
# Convert figure size to inches
figSizeXinches = figSizeXcm/2.54
figSizeYinches = figSizeYcm/2.54

# Padding to move axis labels away from the axis
tickPad = 3
tickLength = 4
markerSize = 4
labelPadY = 3
labelPadX = 3

# Padding around overall figure border (as a fraction of font size)
borderPad = 1.2

# Colours for the line plots - Can use rgb or html
colour1 = 'black'
colour2 = '#D01D3E'
col = ['blue','green','red','cyan','magenta','yellow','black']
# Plot Configuration
plt.rcParams['font.family'] = fontChoice
plt.rcParams['axes.linewidth'] = lineWidth
plt.rcParams["figure.figsize"] = (figSizeXinches, figSizeYinches)

   #### 
fig1 = plt.figure(1)
ax = fig1.add_subplot(111)
fig1.tight_layout(pad=borderPad)

###
ax.plot(D, spatialCorr, label ="omnidirectional", markersize=markerSize,linewidth=lineWidth+0.2,color = 'red', alpha=0.6 )
ax.plot(D, mostCommonCorr, label ="isotropic", markersize=markerSize,linewidth=lineWidth+0.2, color = 'blue' )
ax.set_ylabel(r'$|\rho|$', fontsize=fontSize, labelpad=labelPadY)
ax.set_xlabel(r'$d/\lambda$', fontsize=fontSize, labelpad=labelPadY)
plt.legend(loc='upper right', fontsize = fontSize)



#plt.xlim([0, 12])
plt.xticks(fontsize=fontSize)
#plt.ylim([0, 1.4])
plt.yticks(fontsize=fontSize)
plt.minorticks_off()
ax.tick_params(which='both', direction='in', length=tickLength,
           width=lineWidth,  pad=tickPad, color=colour1)
#ax.grid(color=colour1, linewidth=lineWidth-0.5, alpha=0.1)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')

fig1.savefig('theoretSpatialBasic.eps', format='eps', dpi=200)
    
#TO DO NEXT:
#CHECK https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1330410
#AND https://www.mendeley.com/reference-manager/reader/b6cae3f3-f127-3c5a-b950-2292fefde4da/0fa75d5f-d80b-67ed-29fd-a8f2e8925847
#AND: Good one:https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1043864
#I would say to validate proximity once the correlation is bigger than 0.35 or so. 
#https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1043864
#https://asa.scitation.org/doi/pdf/10.1121/1.1908122
#https://www.mendeley.com/reference-manager/reader/3ee83584-7d05-3989-bd41-8cbc7298998e/3bcf44e2-b152-6641-361c-a44b22faa433
#^^^^^^ FANTASTIC ^^^^^^^^