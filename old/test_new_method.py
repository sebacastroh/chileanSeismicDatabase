import os
import pickle
import random
import numpy as np
import scipy.io as spio
import scipy.integrate as spin
import correctRecords
import matplotlib.pyplot as plt
import multiprocessing
from sklearn.linear_model import LinearRegression

nTest = 10
plt.close('all')
if __name__ == '__main__':
    path = 'events_mat_uncorrected'
    
    with open('p_waves.pkl', 'rb') as f:
        p_waves = pickle.load(f)
    
    n = len(p_waves)
    selected = []
    while len(selected) < nTest:
        i = random.randint(0, n)
        while p_waves[i][4] != 'Valid' or float(p_waves[i][0].split('_')[1][:-2]) < 6.0:
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
        original_station, p_wave = inputs[i]
        new_station = results[i]

        
        g = 9.81
        dt = original_station.dt
        n = len(original_station.acc_1)
        t = np.linspace(0., (n-1)*dt, n)
        
        # fig1 = plt.figure(figsize=(18.76,   6.26))

        # for i in range(3):
        #     acc = getattr(original_station, 'acc_%i' %(i+1))
        #     vel = spin.cumtrapz(acc, dx=dt, initial=0.)
        #     vel -= vel[:p_wave].mean()
        #     dis = spin.cumtrapz(vel, dx=dt, initial=0.)
        #     dis -= dis[:p_wave].mean()

        #     ax = fig1.add_subplot(3,3,i+1)
        #     ax.plot(t, acc/g)
        #     ax.set_ylabel('Acceleration [g]')

        #     ax = fig1.add_subplot(3,3,i+4)
        #     ax.plot(t, vel)
        #     ax.set_ylabel('Velocity [m/s]')

        #     ax = fig1.add_subplot(3,3,i+7)
        #     ax.plot(t, dis*100)
        #     ax.set_ylabel('Displacement [cm]')
        
        # plt.suptitle('Original')
        
        fig2 = plt.figure(figsize=(18.76,   6.26))

        for i in range(3):
            acc = new_station['acc_%i' %(i+1)]
            vel = spin.cumtrapz(acc, dx=dt, initial=0.)
            vel -= vel[:p_wave].mean()
            dis = spin.cumtrapz(vel, dx=dt, initial=0.)
            dis -= dis[:p_wave].mean()
            
            # model = LinearRegression()
            # model.fit(t.reshape(-1,1), dis)
            # trend = model.predict(t.reshape(-1,1))
            
            # dis = dis - trend
            # vel = np.gradient(dis, dt, edge_order=2)
            # acc = np.gradient(vel, dt, edge_order=2)

            ax = fig2.add_subplot(3,3,i+1)
            ax.plot(t, acc/g)
            ax.set_ylabel('Acceleration [g]')

            ax = fig2.add_subplot(3,3,i+4)
            ax.plot(t, vel)
            ax.set_ylabel('Velocity [m/s]')

            ax = fig2.add_subplot(3,3,i+7)
            ax.plot(t, dis*100)
            ax.set_ylabel('Displacement [cm]')

        plt.suptitle('Corrected')

        fig2 = plt.figure(figsize=(18.76,   6.26))

        for i in range(3):
            acc = new_station['acc_%i' %(i+1)]
            vel = spin.cumtrapz(acc, dx=dt, initial=0.)
            vel -= vel[:p_wave].mean()
            dis = spin.cumtrapz(vel, dx=dt, initial=0.)
            dis -= dis[:p_wave].mean()
            
            model = LinearRegression()
            model.fit(t.reshape(-1,1), dis)
            trend = model.predict(t.reshape(-1,1))
            
            dis = dis - trend
            vel = np.gradient(dis, dt, edge_order=2)
            acc = np.gradient(vel, dt, edge_order=2)

            ax = fig2.add_subplot(3,3,i+1)
            ax.plot(t, acc/g)
            ax.set_ylabel('Acceleration [g]')

            ax = fig2.add_subplot(3,3,i+4)
            ax.plot(t, vel)
            ax.set_ylabel('Velocity [m/s]')

            ax = fig2.add_subplot(3,3,i+7)
            ax.plot(t, dis*100)
            ax.set_ylabel('Displacement [cm]')

        plt.suptitle('Corrected')

        plt.draw()
        plt.show()
        plt.pause(5)
        # fig1.savefig(original_station.event_name + '_sta_' + original_station.name + '_original.pdf', bbox_inches='tight', pad_inches=0)
        # fig2.savefig(original_station.event_name + '_sta_' + original_station.name + '_corrected.pdf', bbox_inches='tight', pad_inches=0)
        input("Press Enter to continue...")
        plt.close('all')