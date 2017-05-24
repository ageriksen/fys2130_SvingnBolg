"""this program aims to measure total energy throughout motion of wave"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as man

N = 200 # number of point massses
m = 0.02 #kg mass of single points
k = 10 #kg s^-2 spring constant
mi = np.zeros(N) + m # mass array for mass points
mi[0] = 10000*m
mi[-1] = 10000*m
ki = np.zeros(N - 1) + k # spring constant array for springs
dt = 0.99 * np.sqrt(float(m/k)) # setting timestep. To ensure proper resolution, scaled with angular frequency.
nt = 1200 # number of timesteps
y = np.zeros((N, nt)) # initial position state. Later to be current. box-arry to store time states
i = np.linspace(0, N-1, N) # index array. usable as position along x-axis
y[:,0] = np.sin(7 * np.pi * (i / (N - 1)))
y[:,1] = y[:,0]


for t in range(1, nt-1):
    # time evolution of string
    y[0, t+1] = (dt**2/mi[0])*(-ki[0]*y[0, t] + ki[0]*y[1, t]) + 2*y[0, t] - y[0, t-1]
    y[-1, t+1] =(dt**2/mi[-1])*(-ki[-2]*y[-1, t] + ki[-2]*y[-2, t]) + 2*y[-1, t] - y[-1, t-1]
    for j in range(1, N - 1):
        y[j, t+1] = (dt**2/mi[j])*(-(ki[j-1] + ki[j])*y[j, t] + ki[j-1]*y[j-1, t] + ki[j]*y[j+1, t]) + 2*y[j, t] - y[j, t-1]

#empty energy arrays
Ki = np.zeros((N, nt-1))
Vi = np.zeros((N, nt-1))
E = np.zeros(nt-1)
error = np.zeros(nt-1)

###Calculating energies
# first energies
for i in range(N - 1):
    # Vi[i, 0] = potential(i, 0)
    Vi[i, 0] = 0.5*ki[i]*(y[i+1,0] - y[i, 0])**2
E[0] = np.sum(Vi[:,0])
#last energies
for i in range(N - 1):
    if i < N-1:
        # Vi[i, -1] = potential(i, -2)
        Vi[i, -1] = 0.5*ki[i]*(y[i+1, -2] - y[i, -2])**2
        # Ki[i, -1] = kinetic(i, -2)
        Ki[i, -1] = 0.5*mi[i]*( (y[i, t+1] - y[i, t-1]) / (2*dt) )**2
    else:
        # Ki[i, -1] = kinetic(i, -2)
        Ki[i, -1] = 0.5*mi[i]*( (y[i, t+1] - y[i, t-1]) / (2*dt) )**2
E[-1] = np.sum(Vi[:,-1]) + np.sum(Ki[:,-1])
error[-1] = abs(E[-1] - E[0])/E[0]
# remaining timestate energies
for t in range(1, nt-2):
    for i in range(N):
        if i < N - 1:
            # Vi[i, t] = potential(i, t)
            Vi[i, t] = 0.5*ki[i]*( y[i+1, t] - y[i, t] )**2
            # Ki[i, t] = kinetic(i, t)
            Ki[i, t] = 0.5*mi[i]*( (y[i, t+1] - y[i, t-1]) / (2*dt) )**2
        else:
            # Ki[i, t] = kinetic(i, t)
            Ki[i, t] = 0.5*mi[i]*( (y[i, t+1] - y[i, t-1]) / (2*dt) )**2
    E[t] = np.sum(Vi[:,t]) + np.sum(Ki[:,t])
    error[t] = (abs(E[t]-E[0])) / E[0]

print 'with a near initial energy of %f and a final energy of %f, with an error of %f' %(E[0], E[-1], error[-1])
"""
>> with a near initial energy of 6.069330 and a final energy of 6.075566, with an error of 0.001028
"""

plt.figure()
plt.plot(E[0:-2])
plt.ylim(0, np.max(E)+ 0.2)
plt.xlabel('timesteps i')
plt.ylabel('energy of system')
plt.title('Total energy of system over time')

plt.figure()
plt.plot(error[0:-2])
plt.ylim(-1, 1)
plt.xlabel('timesteps i')
plt.ylabel('relative difference scaled to initial energy')
plt.title('relative error in energy loss for system over time')

plt.show()