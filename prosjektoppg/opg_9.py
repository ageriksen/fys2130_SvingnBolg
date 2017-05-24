"""program to expand on previous iterations of simulated waves on a string"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as man

N = 400
m = 0.02 #kg mass of single points
k = 10 #kg s^-2 spring constant
dt =  1 * np.sqrt(float(m/k)) # setting timestep. To ensure proper resolution, scaled with angular frequency.
nt = 1200 # number of timesteps
##setting variables in arrays
mi = np.zeros(N) + m # mass array for mass points
mi[0] = 10000*m
for j in range(N-1):
    if 200 <= j:
        mi[j] = 3*m
mi[-1] = 30000*m
ki = np.zeros(N - 1) + k # spring constant array for springs
y = np.zeros((N, nt)) # initial position state. Later to be current. box-arry to store time states
i = np.linspace(0, N-1, N) # index array. usable as position along x-axis

# setting initial wave (currently making the straight triangle)
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

# running through timedevelopement
for t in range(1, nt-1):
    y[0, t+1] = (dt**2/mi[0])*(-ki[0]*y[0, t] + ki[0]*y[1, t]) + 2*y[0,t] - y[0, t-1]
    y[-1, t+1] = (dt**2/mi[-1])*(-ki[-2]*y[-1, t] + ki[-2]*y[-2, t]) + 2*y[-1, t] - y[-1, t-1]
    for j in range(1, N-2):
        y[j, t+1] = (dt**2/mi[j])*(-(ki[j-1] + ki[j])*y[j,t] + ki[j-1]*y[j-1,t] + ki[j]*y[j+1,t]) + 2*y[j,t] - y[j,t-1]

## read off of plot. maxima and minima for incomin, transmitted and reflected wave.
#y[np.argmax(y[140:170, t]), t]
Ai = np.max(y[126:132, 100])
Ar = np.min(y[76:85, 290])
At = np.max(y[264:271, 290])
print 'Ai ', Ai
print 'Ar ', Ar
print 'At ', At
#from modifying formula for accoustinc impedance
Z0 = np.sqrt(mi[149]*ki[149])
Z1 = np.sqrt(mi[249]*ki[249])
Z1O2_num = At/(Ai - Ar)
Z1O2_ana = Z0/Z1
ZError = (abs(Z1O2_num - Z1O2_ana))/Z1O2_ana
print 'The measured relation of impedance is %.3f as opposed to the analytic %.3f with an error %.3f' %(Z1O2_num, Z1O2_ana, ZError)
"""
>> The measured relation of impedance is 0.555 as opposed to the analytic 0.577 with an error 0.038
"""

#saving array
# np.save('TriWaveLongrope.npy', y)

#ploting single timesstates
plt.figure()
plt.plot(y[:,100])
plt.figure()
plt.plot(y[:,290])
plt.show()

###animating
#setup
fig = plt.figure()
ax = plt.axes(xlim=(-1, N+1), ylim=(-2, 2))
line, = ax.plot([], [])
#background
def init():
    line.set_data([], [])
    return line,
#animated func
def animate(j):
    line.set_data(i, y[:,j])
    return line,
#animator
anim = man.FuncAnimation(fig, animate, init_func=init, frames=1200, interval=20, blit=True)
plt.show()