import os
import pickle
import numpy as np
import scipy.integrate as spin

path = os.path.abspath('./events_pkl_corrected_v2/')
filenames = os.listdir(path)

d1 = []
d2 = []

for filename in filenames:
    full_path = os.path.join(path, filename)
    with open(full_path, 'rb') as f:
        event = pickle.load(f)

    for station in event.values():
        dt = station['dt']
        for i in range(3):
            key = 'acc_%i' %(i+1)
            acc = station[key]
            if np.allclose(acc, np.zeros_like(acc)):
                continue
            ia = spin.cumtrapz(acc**2, dx=dt, initial=0.)
            ia /= ia[-1]
            p1 = np.argmin(np.abs(ia - 0.05))
            p2 = np.argmin(np.abs(ia - 0.75))
            p3 = np.argmin(np.abs(ia - 0.95))

            d1.append((p2-p1)*dt)
            d2.append((p3-p1)*dt)

np.save('d5_75.npy', np.array(d1))
np.save('d5_95.npy', np.array(d2))
