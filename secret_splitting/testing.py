import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
 
#Simulation window parameters
xMin=-2;xMax=2;
yMin=-2;yMax=2;
xDelta=xMax-xMin;yDelta=yMax-yMin; #rectangle dimensions
areaTotal=xDelta*yDelta;
 
#Point process parameters
lambda0 = .4; #intensity (ie mean density) of the Poisson process
 
#Simulate Poisson point process
numbPoints = scipy.stats.poisson( lambda0*areaTotal ).rvs()#Poisson number of points
xx = xDelta*scipy.stats.uniform.rvs(0,1,((numbPoints,1)))+xMin#x coordinates of Poisson points
yy = yDelta*scipy.stats.uniform.rvs(0,1,((numbPoints,1)))+yMin#y coordinates of Poisson points
#Plotting
fignum = 1
fig = plt.figure(fignum, figsize=(8, 5))
plt.scatter(xx,yy, edgecolor='b', facecolor='none', alpha=0.5 )
plt.xlabel("x"); plt.ylabel("y")
plt.text(1, 0,"A1", fontsize=12)
plt.scatter([0, -1, 1], [0, 0, 0], s = 5000, c=(0,0,0), alpha=0.1)
plt.text(-1, 0,"A2", fontsize=12)
plt.text(0, 0,"B", fontsize=12)
plt.show()
plt.close()

breakpoint()