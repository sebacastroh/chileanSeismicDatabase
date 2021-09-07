import sys
import numpy as np
import scipy.io as spio
import scipy.integrate as spin
import matplotlib.pyplot as plt

event_name = sys.argv[1]

event = spio.loadmat(event_name, struct_as_record=False, squeeze_me=True)

event.pop('__version__')
event.pop('__header__')
event.pop('__globals__')

for station in event.itervalues():
    dt = station.dt
    n = len(station.acc_1)
    t = np.linspace(0., (n-1)*dt, n)
    
    plt.figure()
    plt.subplot(311)
    plt.plot(t, station.acc_1)

    vel = spin.cumtrapz(station.acc_1, dx=dt, initial=0.)
    plt.subplot(312)
    plt.plot(t, vel)

    dis = spin.cumtrapz(vel, dx=dt, initial=0.)
    plt.subplot(313)
    plt.plot(t, dis)

    plt.figure()
    plt.subplot(311)
    plt.plot(t, station.acc_2)

    vel = spin.cumtrapz(station.acc_2, dx=dt, initial=0.)
    plt.subplot(312)
    plt.plot(t, vel)

    dis = spin.cumtrapz(vel, dx=dt, initial=0.)
    plt.subplot(313)
    plt.plot(t, dis)

    plt.figure()
    plt.subplot(311)
    plt.plot(t, station.acc_3)

    vel = spin.cumtrapz(station.acc_3, dx=dt, initial=0.)
    plt.subplot(312)
    plt.plot(t, vel)

    dis = spin.cumtrapz(vel, dx=dt, initial=0.)
    plt.subplot(313)
    plt.plot(t, dis)



plt.ion()
plt.show()
