import os
import pickle
import random
import numpy as np
import scipy.io as spio
import scipy.integrate as spin
import correctRecords
import matplotlib.pyplot as plt
import multiprocessing

nTest = 10

if __name__ == '__main__':
    path = 'events_mat_uncorrected'
    
    with open('p_waves.pkl', 'rb') as f:
        p_waves = pickle.load(f)
    
    n = len(p_waves)
    selected = []
    while len(selected) < nTest:
        i = random.randint(0, n)
        while p_waves[i][4] != 'Valid':
            i = random.randint(0, n)
    
        event = p_waves[i][0]
        station = p_waves[i][1]
        p_loc = p_waves[i][3]
        group = [event, station, p_loc]
        
        if group not in selected:
            selected.append(group)
            
    pool = multiprocessing.Pool(processes=nTest)
    
    inputs = []
    for s in selected:
        event, station, p_loc = s
        mat_event = spio.loadmat(os.path.join(path, event), struct_as_record=False, squeeze_me=True)
        mat_station = mat_event[station]
        inputs.append([mat_station, p_loc])
    
    results = pool.starmap(correctRecords.correctStation, inputs)
    pool.close()

    for i in range(nTest):
        original_station, p_loc = inputs[i]
        new_station = results[i]

        t = np.linspace(0., (len(original_station.acc_1)-1)*mat_station.dt, len(original_station.acc_1))

        fig1 = plt.figure()

        a1a = fig1.add_subplot(331)
        a1a.plot(t, original_station.acc_1)

        a2a = fig1.add_subplot(332)
        a2a.plot(t, original_station.acc_2)

        a3a = fig1.add_subplot(333)
        a3a.plot(t, original_station.acc_3)

        a1v = fig1.add_subplot(334)
        a1v.plot(t, spin.cumtrapz(original_station.acc_1, dx=original_station.dt, initial=0.))

        a2v = fig1.add_subplot(335)
        a2v.plot(t, spin.cumtrapz(original_station.acc_2, dx=original_station.dt, initial=0.))

        a3v = fig1.add_subplot(336)
        a3v.plot(t, spin.cumtrapz(original_station.acc_3, dx=original_station.dt, initial=0.))

        a1d = fig1.add_subplot(337)
        a1d.plot(t, spin.cumtrapz(spin.cumtrapz(original_station.acc_1, dx=original_station.dt, initial=0.), dx=original_station.dt, initial=0.))

        a2d = fig1.add_subplot(338)
        a2d.plot(t, spin.cumtrapz(spin.cumtrapz(original_station.acc_2, dx=original_station.dt, initial=0.), dx=original_station.dt, initial=0.))

        a3d = fig1.add_subplot(339)
        a3d.plot(t, spin.cumtrapz(spin.cumtrapz(original_station.acc_3, dx=original_station.dt, initial=0.), dx=original_station.dt, initial=0.))

        fig2 = plt.figure()

        a1a = fig2.add_subplot(331)
        a1a.plot(t, new_station['acc_1'])

        a2a = fig2.add_subplot(332)
        a2a.plot(t, new_station['acc_2'])

        a3a = fig2.add_subplot(333)
        a3a.plot(t, new_station['acc_3'])

        a1v = fig2.add_subplot(334)
        a1v.plot(t, spin.cumtrapz(new_station['acc_1'], dx=new_station['dt'], initial=0.))

        a2v = fig2.add_subplot(335)
        a2v.plot(t, spin.cumtrapz(new_station['acc_2'], dx=new_station['dt'], initial=0.))

        a3v = fig2.add_subplot(336)
        a3v.plot(t, spin.cumtrapz(new_station['acc_3'], dx=new_station['dt'], initial=0.))

        a1d = fig2.add_subplot(337)
        a1d.plot(t, spin.cumtrapz(spin.cumtrapz(new_station['acc_1'], dx=new_station['dt'], initial=0.), dx=new_station['dt'], initial=0.))

        a2d = fig2.add_subplot(338)
        a2d.plot(t, spin.cumtrapz(spin.cumtrapz(new_station['acc_2'], dx=new_station['dt'], initial=0.), dx=new_station['dt'], initial=0.))

        a3d = fig2.add_subplot(339)
        a3d.plot(t, spin.cumtrapz(spin.cumtrapz(new_station['acc_3'], dx=new_station['dt'], initial=0.), dx=new_station['dt'], initial=0.))

        plt.draw()
        plt.show()
        plt.pause(5)
        input("Press Enter to continue...")
        plt.close('all')