import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as man

i = np.linspace(0, N-1, N) # index array. usable as position along x-axis

y_stand = np.load('timestates.npy')

### animating
fig = plt.figure()
ax = plt.axes(xlim=(-0.5, 200.5), ylim=(-1.2, 1.2))
line, = ax.plot([], [])
#init func
def init():
    line.set_data([],[])
    return line,

#animating function
def animate(j):
    line.set_data(i, y_stand[:, j])
    return line,

anim = man.FuncAnimation(fig, animate, init_func=init, frames=200, interval=30, blit=True)
plt.show()
