"""This program aims to simulate oscillations in a string"""

from math import pi, sin
import matplotlib.pyplot as plt
import numpy as np

#initializing variables in arrays
N = 200 # number of timesteps
m = 0.02 #kg mass of single points
k = 10 #kg s^-2 spring constant
dt = 0.8 * np.sqrt(float(m/k)) # setting timestep. To ensure proper resolution, scaled with angular frequency.
nt = 1200 # number of timesteps
mi = np.zeros(N) + m # mass array for mass points
ki = np.zeros(N - 1) + k # spring constant array for springs
y0 = np.zeros((N, nt)) # initial position state. Later to be current. box-arry to store time states
i = np.linspace(0, N-1, N) # index array. usable as position along x-axis
y0[:,0] = np.sin(7 * pi * (i / (N - 1)))
print i
# setting initial states. static start, ensures previous state equal.
y_l = y0[:,0] # "y-last"


# initializing "current" array
y_n = np.zeros(N) # "y-next" next position state

#nested loop updating all points in y per time
# for t in range(1):
#     y_n0
#     y_nN
#     for i in range(1, N-1):
#         y_n[i]
#     y_l = y_0
#     y_0 = y_n
#
# plt.plot(i, y_0)