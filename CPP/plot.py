#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

data1 = np.loadtxt('u.txt')

m = int(data1[0])   # number of nodes
t = int(data1[1])   # number of time steps

data2 = data1[2:]

u_sol = data2.reshape(t, m)

for i in xrange(0, t):
    plt.plot(u_sol[i])


plt.show()
