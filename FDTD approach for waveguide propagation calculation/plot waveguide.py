# -*- coding: utf-8 -*-
"""
Created on Sun Sep  9 22:35:18 2018

@author: zdlls
"""

import matplotlib.pyplot as plt
import numpy as np

data=np.loadtxt('waveguide.txt')
x=data[:,0]
y=data[:,1]
z=data[:,2]
plt.tripcolor(x,y,z)
plt.show()