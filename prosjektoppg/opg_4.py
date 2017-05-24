"""This program aims to simulate oscillations in a string"""
import matplotlib.pyplot as plt
import matplotlib.animation as man
import numpy as np

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
y[:,0] = np.sin(7 * np.pi * (i / (N - 1)))
y[:,1] = y[:,0]


for t in range(1, nt-1):
    y[0, t+1] = (dt**2/mi[0])*(-ki[0]*y[0, t] + ki[0]*y[1, t]) + 2*y[0, t] - y[0, t-1]
    y[-1, t+1] =(dt**2/mi[-1])*(-ki[-2]*y[-1, t] + ki[-2]*y[-2, t]) + 2*y[-1, t] - y[-1, t-1]
    for j in range(1, N - 1):
        y[j, t+1] = (dt**2/mi[j])*(-(ki[j-1] + ki[j])*y[j, t] + ki[j-1]*y[j-1, t] + ki[j]*y[j+1, t]) + 2*y[j, t] - y[j, t-1]

###saving y as npy file
np.save('timestates.npy', y)

### plotting animation of waves
#setup fig, axes and plot element to animate
fig = plt.figure()
ax = plt.axes(xlim=(-1, N+1), ylim=(-2, 2))
line, = ax.plot([], [])
#init func, plot background
def init():
    line.set_data([], [])
    return line,
#sequentially called animation function
def animate(j):
    line.set_data(i, y[:,j])
    return line,
#call of animator "blit=True" only redraws change.
anim = man.FuncAnimation(fig, animate, init_func=init, frames=1200, interval=20, blit=True)
plt.show()

# ###ploting exerpts
plt.figure()
plt.plot(i, y[:,0])
plt.xlabel('position i')
plt.ylabel('position y')
plt.title('state at timestep nt[i]=0')
plt.ylim(-1.2, 1.2)
plt.figure()
plt.plot(i, y[:,15])
plt.xlabel('position i')
plt.ylabel('position y')
plt.title('state at timestep nt[i]=15')
plt.ylim(-1.2, 1.2)
plt.figure()
plt.plot(i, y[:, 30])
plt.xlabel('position i')
plt.ylabel('position y')
plt.title('state at timestep nt[i]=30')
plt.ylim(-1.2, 1.2)
plt.show()