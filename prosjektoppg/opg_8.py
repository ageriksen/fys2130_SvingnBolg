""" this program aims to simulate the propagation of a single exitation on a string. 
    centered on point 30 to start off. We want to send it all going to the right. 
    This means g(x,t) = G(xt - wt). """

# from opg_4 import i, y, nt, mi, ki, dt, N
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as man

N = 200 # number of point massses
m = 0.02 #kg mass of single points
k = 10 #kg s^-2 spring constant
dt = 0.99 * np.sqrt(float(m/k)) # setting timestep. To ensure proper resolution, scaled with angular frequency.
nt = 1200 # number of timesteps
ki = np.zeros(N - 1) + k # spring constant array for springs

mi = np.zeros(N) + m # mass array for mass points
mi[0] = 10000*m
mi[-1] = 10000*m

i = np.linspace(0, N-1, N) # index array. usable as position along x-axis
y = np.zeros((N, nt)) # initial position state. Later to be current. box-arry to store time states

for j in range(len(i)/2):
    if 0 <= j <= 29:
        y[j, 0] = (j)/30.
    elif 30 <= j <= 58:
        y[j, 0] = (59 - j)/30.

for j in range(len(i)/2):
    if 1 <= j <= 30:
        y[j, 1] = (j-1)/30.
    elif 31 <= j <= 59:
        y[j, 1] = (60 - j)/30.

# y[:, 1] = y[:, 0]

if __name__ == '__main__':
    for t in range(1, nt-1):
        y[0, t+1] = (dt**2/mi[0])*(-ki[0]*y[0, t] + ki[0]*y[1, t]) + 2*y[0, t] - y[0, t - 1]
        y[-1, t+1] = (dt**2/mi[-1])*(-ki[-2]*y[-1, t] + ki[-2]*y[-2, t]) + 2*y[-1, t] - y[-1, t - 1]
        for j in range(1, N-1):
            y[j, t+1] = (dt**2/mi[j])*(-(ki[j - 1] + ki[j])*y[j, t] + ki[j - 1]*y[j - 1, t] + ki[j]*y[j + 1, t]) \
                        + 2*y[j, t] - y[j, t - 1]

    ###save new y as npy file
    # np.save('triangle_wave.npy', y)

    ###plotting animation
    #setup fig, axes and plot-element to animate
    fig = plt.figure()
    ax = plt.axes(xlim=(-.5, N+.5), ylim=(-1.5, 1.5))
    line, = ax.plot([], [])
    #init function to plot background
    def init():
        line.set_data([],[])
        return line,
    #sequentially called animate func
    def animate(j):
        line.set_data(i, y[:,j])
        return line,

    #animator call, blit=true only draw changes.
    anim = man.FuncAnimation(fig, animate, init_func=init, frames=1200, interval=20, blit=True)
    plt.show()

    plt.figure()
    plt.plot(y[:, 50])
    plt.xlabel('position in distance, i [dimless]')
    plt.ylabel('position in exitation y [dimless]')
    plt.title('Triangle wave propagating to the right, [t=50]')
    plt.show()