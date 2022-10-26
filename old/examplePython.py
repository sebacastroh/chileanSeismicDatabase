import pickle
import numpy as np
import scipy.integrate as spin
import matplotlib.pyplot as plt

EarthquakeName = '20171013_4Mw_20.29S_69.14W_91KM.pkl'

with open(EarthquakeName, 'rb') as f:
    event = pickle.load(f)

print('##################################')
print('Event: ' + EarthquakeName)
print('Number of stations: %i' %len(event))
print('##################################')

# Let's plot the acceleration, velocity and displacement of the first component for the first station
station = event['st00']

print('Attributes station "st00":')
for key,value in station.items():
    print('station["' + key + '"] = ', value)

fig = plt.figure()
a1 = fig.add_subplot(311)

t = np.linspace(0., (len(station['acc_1'])-1)*station['dt'], len(station['acc_1']))

a1.plot(t, station['acc_1']/9.81)
a1.set_ylabel('a(t) [g]')

a2 = fig.add_subplot(312)
v = spin.cumtrapz(station['acc_1'], dx=station['dt'], initial=0.)
a2.plot(t, v)
a2.set_ylabel('v(t) [m/s]')

a3 = fig.add_subplot(313)
d = spin.cumtrapz(v, dx=station['dt'], initial=0.)
a3.plot(t, d)
a3.set_ylabel('d(t) [m]')
a3.set_xlabel('Time [s]')

fig.suptitle('Station ' + station['name'] + '\nComponent ' + station['component_1'])
plt.show()
