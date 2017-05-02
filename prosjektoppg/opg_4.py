"""This program aims to simulate oscillations in a string"""

import matplotlib.pyplot as plt
import matplotlib.animation as man
import numpy as np

##defining functions

def Extremal(i, t):
    if i == 0:
        return (dt**2/mi[i])*(-ki[i]*y[i, t] + ki[i]*y[i+1, t]) + 2*y[i, t] - y[i, t-1]
    elif i == N-1:
        return (dt**2/mi[i])*(-ki[i-1]*y[i, t] + ki[i-1]*y[i-1, t]) + 2*y[i, t] - y[i, t-1]
    else:
        print "input wrong position"
        quit()

def NextY(i, t):
    return (dt**2/mi[i])*( -(ki[i-1]+ki[i])*y[i, t] + ki[i-1]*y[i-1, t] + ki[i]*y[i+1, t]) + 2*y[i, t] - y[i, t-1]

##initializing variables
N = 200 # number of timesteps
m = 0.02 #kg mass of single points
k = 10 #kg s^-2 spring constant
dt = 0.8 * np.sqrt(float(m/k)) # setting timestep. To ensure proper resolution, scaled with angular frequency.
nt = 1200 # number of timesteps
##setting variables in arrays
mi = np.zeros(N) + m # mass array for mass points
ki = np.zeros(N - 1) + k # spring constant array for springs
y = np.zeros((N, nt)) # initial position state. Later to be current. box-arry to store time states
i = np.linspace(0, N-1, N) # index array. usable as position along x-axis
y[:,0] = np.sin(7 * np.pi * (i / (N - 1)))

### manual 1st step
y[0, 1] = Extremal(0, 1)
y[-1, 1] = Extremal(N-1, 1)
for j in range(1, N-1):
    y[j, 1] = NextY(j, 1)
    
### finishing from 2nd timestep and on
for t in range(1, nt-1):
    y[0, t] = Extremal(0, t)
    y[-1, t] = Extremal(N - 1, t)
    for j in range(1, N - 1):
        y[j, t] = NextY(j, t)

### plotting animation of waves