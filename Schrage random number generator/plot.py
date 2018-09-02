# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import pylab as pl

data=np.loadtxt('D01.txt')
pl.plot(data[:,0],data[:,1],'r.')
pl.xlabel('X Data')
pl.ylabel('Y Data')
pl.xlim(-0.2,1.2)
pl.ylim(-0.2,1.2)
pl.show()
