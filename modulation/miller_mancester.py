# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 09:41:46 2020

@author: chrys
"""


import numpy as np
import matplotlib.pyplot as plt
import encoding as en

binary_message = [1,0,0,1,1,0,0,1,0,1,1,0,1,0,1]

bb = np.repeat(binary_message,2)
bb = np.append(bb,binary_message[-1])
manch = en.manchester(binary_message)
miller = en.miller(binary_message)
fm0 = en.fm0(binary_message)


plt.subplot(4,1,1)
plt.title('Message Signal')
plt.plot(bb, marker='', color='black', drawstyle='steps-post')
#plt.ylabel('Amplitude')
#plt.xlabel('t')
plt.xticks([], [])
plt.yticks([],[])

             
for i in range(len(binary_message)):
    plt.text(2*i+1, 0.5 , str(binary_message[i]), fontsize=14)

plt.subplot(4,1,2)
plt.title('fm0 code')
plt.plot(fm0, marker='', color='red', drawstyle='steps-post')
#plt.ylabel('Amplitude')
#plt.xlabel('t')
plt.xticks([], [])
plt.yticks([],[])

plt.subplot(4,1,3)
plt.title('Manchester code')
plt.plot(manch, marker='', color='blue', drawstyle='steps-post')
#plt.ylabel('Amplitude')
#plt.xlabel('t')
plt.xticks([], [])
plt.yticks([],[])

plt.subplot(4,1,4)
plt.title('Miller code')
plt.plot(miller, marker='', color="purple", drawstyle='steps-post')
#plt.ylabel('Amplitude')
#plt.xlabel('t')
plt.xticks([], [])
plt.yticks([],[])

print(binary_message)

plt.subplots_adjust(hspace=1)
plt.rc('font', size=15)
fig = plt.gcf()
fig.set_size_inches(16, 9)

#fig.savefig('Amplitude Modulation.png', dpi=100)