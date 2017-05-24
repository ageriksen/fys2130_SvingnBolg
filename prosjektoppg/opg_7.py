""" this program aims to simulate the propagation of a single exitation on a string. """

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as man

##initializing variables
N = 200 # number of point massses
m = 0.02 #kg mass of single points
k = 10 #kg s^-2 spring constant
dt = 0.99 * np.sqrt(float(m/k)) # setting timestep. To ensure proper resolution, scaled with angular frequency.
nt = 1200 # number of timesteps
##setting variables in arrays
mi = np.zeros(N) + m # mass array for mass points
mi[0] = 10000*m
mi[-1] = 10000*m
ki = np.zeros(N - 1) + k # spring constant array for springs
y = np.zeros((N, nt)) # initial position state. Later to be current. box-arry to store time states
i = np.linspace(0, N-1, N) # index array. usable as position along x-axis

for j in range(len(i)):
    if 70 <= j <= 99:
        y[j, 0] = (j - 69)/30.
    elif 100 <= j <= 128:
        y[j, 0] = (129 - j)/30.
    else:
        y[j, 0] = 0
y[:, 1] = y[:, 0]

## running time developement
for t in range(1, nt-1):
    y[0, t+1] = (dt**2/mi[0])*(-ki[0]*y[0, t] + ki[0]*y[1, t]) + 2*y[0, t] - y[0, t - 1]
    y[-1, t+1] = (dt**2/mi[-1])*(-ki[-2]*y[-1, t] + ki[-2]*y[-2, t]) + 2*y[-1, t] - y[-1, t - 1]
    for j in range(1, N-1):
        y[j, t+1] = (dt**2/mi[j])*(-(ki[j - 1] + ki[j])*y[j, t] + ki[j - 1]*y[j - 1, t] + ki[j]*y[j + 1, t]) \
                    + 2*y[j, t] - y[j, t - 1]

###save new y as npy file
# np.save('triangle_wave.npy', y)

### ploting snapshots of motion
plt.figure()
plt.plot(y[:,0])
plt.xlabel('position i')
plt.ylabel('excitation amplitude. [dimensionless]')
plt.title('Triangle wave at time [t=0]')

plt.figure()
plt.plot(y[:,50])
plt.ylim(-.1, 1)
plt.xlabel('position i')
plt.ylabel('excitation amplitude. [dimensionless]')
plt.title('Triangle wave at time [t=50]')
plt.show()

###plotting animation
#setup fig, axes and plot-element to animate
# fig = plt.figure()
# ax = plt.axes(xlim=(-.5, N+.5), ylim=(-1.5, 1.5))
# line, = ax.plot([], [])
# #init function to plot background
# def init():
#     line.set_data([],[])
#     return line,
# #sequentially called animate func
# def animate(j):
#     line.set_data(i, y[:,j])
#     return line,
#
# #animator call, blit=true only draw changes.
# anim = man.FuncAnimation(fig, animate, init_func=init, frames=1200, interval=20, blit=True)
# plt.show()
# # vi ser bolgen splitter seg, og beveger seg i hver sin retning. Den reflekteres av de faste endepunktene, og
# #gar sammen igjen mot slutten av simulasjonnen. Dette vil repeteres.