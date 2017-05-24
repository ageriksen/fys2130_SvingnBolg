"""this program aims to plot 3 periods of the hundredth element of the standing wave on the string."""

import numpy as np
import matplotlib.pyplot as plt

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

###calculating numerical period, comparing.
f_expected = np.sqrt(k/m)*(3.5/200) #expected frequency from plots of prob. 4
y_ma1 = 0 # variables to store maxima
y_ma2 = 0
y_m1 = 0 #variables to store indices of maxima
y_m2 = 0 # ^^^

t_10 = 10/f_expected # time passed over 10 periods.
nT_10 = int(t_10/dt)+1 # picking 1 more timestep than the 10 periods I want
y99 = y[99, 0:nT_10] # picking out the 10 periods from our expected frequency

# plotting
plt.plot(y99)
plt.xlabel('timesteps')
plt.ylabel('vertical displacement')
plt.title('vertical displacement 100th point over time')
plt.show()

for i in range(len(y99)):
    if 27 < i < 30: # these limits read from plot
        if y99[i] > y_ma1:
            y_ma1 = y99[i]
            y_m1 = i
    elif 85 < i < 88:
        if y99[i] > y_ma2:
            y_ma2 = y99[i]
            y_m2 = i
    else:
        pass
T = (y_m2 - y_m1)*dt
f = 1/T
error = (abs(f-f_expected))/f_expected
print 'we have an error of %.3f' %(error), ' from our expected freq. %.3f' %(f_expected), ' and our evaluated %.3f' %(f)
"""
>> we have an error of 0.005  from our expected freq. 0.391  and our evaluated 0.389
"""

# FFT of y99
FFTy99 = np.abs(np.fft.fft(y99))**2
# freq = np.fft.fftfreq(len(y99), dt)
fr = np.linspace(0, (1/(2*dt)), len(y99)/2)
fftMax =fr[np.argmax(FFTy99[0:11])] # unnecessary to consider higher frequencies than 1s^-1
fftError = (abs(fftMax - f_expected))/f_expected
print 'the fft gives a spike at %.3f, which means an error of %.3f with regards to expected' %(fftMax, fftError)
"""
>> the fft gives a spike at 0.392, which means an error of 0.002 with regards to expected
"""
plt.plot(fr, FFTy99[0:(len(FFTy99)/2)])
plt.xlabel('frequency [s^-1]')
plt.ylabel('arbitrary units, Fourier Coeff.')
plt.title('fft of frequencies present in oscillations of 100th element')
plt.show()