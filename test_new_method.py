import os
import pickle
import random
import numpy as np
import scipy.io as spio
import correctRecords
import matplotlib.pyplot as plt

path = 'events_mat_uncorrected'

with open('p_waves.pkl', 'rb') as f:
    p_waves = pickle.load(f)


n = len(p_waves)
i = random.randint(0, n)
while p_waves[i][4] != 'Valid':
    i = random.randint(0, n)

event = p_waves[i][0]
station = p_waves[i][1]
p_loc = p_waves[i][3]


mat_event = spio.load(os.path.join(path, event), struct_as_record=False, squeeze_me=True)
mat_station = mat_event[station]

t = np.linspace(0., (len(mat_station.acc_1)-1)*mat_station.dt, len(mat_station.acc_1))

new_station = correctRecords.correctStation(new_station, p_wave)

# fig = plt.figure()
# a1 = fig.add_subplot(311)
# a1.plot(t, mat_station.acc_1)

# a1 = fig.add_subplot(312)
# a1.plot(t, mat_station.acc_1)

# a1 = fig.add_subplot(313)
# a1.plot(t, mat_station.acc_1)