import numpy as np
import scipy.io as spio
import matplotlib.pyplot as plt
from myCodes import SeismicLibrary

codes = spio.loadmat('spectra.mat', struct_as_record=False, squeeze_me=True)
event = spio.loadmat('./events_mat_corrected_v2/20100227_8.8Mw_36.10S_73.08W_30KM.mat', struct_as_record=False, squeeze_me=True)

event.pop('__version__')
event.pop('__globals__')
event.pop('__header__')

Tn = np.logspace(-2, 2, 1001)
xi = 0.05
g = 9.81

fig1, (axb1, axc1, axd1) = plt.subplots(3,1)
fig2, (axb2, axc2, axd2) = plt.subplots(3,1)
fig3, (axb3, axc3, axd3) = plt.subplots(3,1)
fig4, (axb4, axc4, axd4) = plt.subplots(3,1)
fig5, (axb5, axc5, axd5) = plt.subplots(3,1)
#fig6, ax = plt.subplots(1,1)

fig1.set_size_inches([6.4, 6.51])
fig2.set_size_inches([6.4, 6.51])
fig3.set_size_inches([6.4, 6.51])
fig4.set_size_inches([6.4, 6.51])
fig5.set_size_inches([6.4, 6.51])

kb = True
kc = True
kd = True

ratio_b = np.empty((0,1001))
ratio_c = np.empty((0,1001))
ratio_d = np.empty((0,1001))

for station in event.values():
    spectra_rot = SeismicLibrary.SpectraRot(station.acc_1, station.acc_2, station.dt, Tn, xi)

    Sax = spectra_rot[0]
    Say = spectra_rot[90]
    Srotd0 = np.min(spectra_rot, axis=0)
    Srotd50 = np.median(spectra_rot, axis=0)
    Srotd100 = np.max(spectra_rot, axis=0)

    if station.Vs30 < 180:
        print('suelo e')
    elif station.Vs30 < 350:
        if kb:
            label = 'Registros Maule 2010'
            kb = False
        else:
            label = ''

        axd1.semilogx(Tn, Srotd50/g, c='gray', alpha=0.5, label=label)
        axd2.semilogx(Tn, Srotd100/g, c='gray', alpha=0.5, label=label)

        axd3.semilogx(Tn,Tn*Tn*Srotd50/(4.*np.pi**2), c='gray', alpha=0.5, label=label)
        axd4.semilogx(Tn, Tn*Tn*Srotd100/(4.*np.pi**2), c='gray', alpha=0.5, label=label)

        axd5.semilogx(Tn, Srotd100/Srotd50, c='gray', alpha=0.5, label=label)

        ratio_d = np.vstack((ratio_d, Srotd100/Srotd50))

    elif station.Vs30 < 500:
        if kc:
            label = 'Registros Maule 2010'
            kc = False
        else:
            label = ''

        axc1.semilogx(Tn, Srotd50/g, c='gray', alpha=0.5, label=label)
        axc2.semilogx(Tn, Srotd100/g, c='gray', alpha=0.5, label=label)

        axc3.semilogx(Tn,Tn*Tn*Srotd50/(4.*np.pi**2), c='gray', alpha=0.5, label=label)
        axc4.semilogx(Tn, Tn*Tn*Srotd100/(4.*np.pi**2), c='gray', alpha=0.5, label=label)

        axc5.semilogx(Tn, Srotd100/Srotd50, c='gray', alpha=0.5, label=label)

        ratio_c = np.vstack((ratio_c, Srotd100/Srotd50))

    elif station.Vs30 < 900:
        if kd:
            label = 'Registros Maule 2010'
            kd = False
        else:
            label = ''

        axb1.semilogx(Tn, Srotd50/g, c='gray', alpha=0.5, label=label)
        axb2.semilogx(Tn, Srotd100/g, c='gray', alpha=0.5, label=label)

        axb3.semilogx(Tn,Tn*Tn*Srotd50/(4.*np.pi**2), c='gray', alpha=0.5, label=label)
        axb4.semilogx(Tn, Tn*Tn*Srotd100/(4.*np.pi**2), c='gray', alpha=0.5, label=label)

        axb5.semilogx(Tn, Srotd100/Srotd50, c='gray', alpha=0.5, label=label)

        ratio_b = np.vstack((ratio_b, Srotd100/Srotd50))

    else:
        print('suelo a')

Tn_433 = codes['T']
Tn_2745 = codes['T']

pos = np.argmin(np.abs(Tn_433-5))

axb1.semilogx(Tn_433, codes['Sa_NCh433'][1,0], label='NCh433+DS61 - Zona 2')
axc1.semilogx(Tn_433, codes['Sa_NCh433'][1,1], label='NCh433+DS61 - Zona 2')
axd1.semilogx(Tn_433, codes['Sa_NCh433'][1,2], label='NCh433+DS61 - Zona 2')

axb1.semilogx(Tn_2745, codes['Sa_NCh2745'][1,0], label='NCh2745 - Zona 2')
axc1.semilogx(Tn_2745, codes['Sa_NCh2745'][1,1], label='NCh2745 - Zona 2')
axd1.semilogx(Tn_2745, codes['Sa_NCh2745'][1,2], label='NCh2745 - Zona 2')

axb1.text(0.05, 0.8, 'Suelo Tipo B', transform = axb1.transAxes)
axc1.text(0.05, 0.8, 'Suelo Tipo C', transform = axc1.transAxes)
axd1.text(0.05, 0.8, 'Suelo Tipo D', transform = axd1.transAxes)

axb1.set_ylabel('Sa [g]')
axc1.set_ylabel('Sa [g]')
axd1.set_ylabel('Sa [g]')

axb1.legend(loc=1)
axc1.legend(loc=1)
axd1.legend(loc=1)

axb1.set_xlim([1e-1, 10])
axc1.set_xlim([1e-1, 10])
axd1.set_xlim([1e-1, 10])

#axb1.set_ylim([0, 5.13])
#axc1.set_ylim([0, 5.13])
#axd1.set_ylim([0, 5.13])

axd1.set_xlabel('Periodo [s]')
fig1.suptitle('Aceleración espectral RotD50')

axb2.semilogx(Tn_433, codes['Sa_NCh433'][1,0], label='NCh433+DS61 - Zona 2')
axc2.semilogx(Tn_433, codes['Sa_NCh433'][1,1], label='NCh433+DS61 - Zona 2')
axd2.semilogx(Tn_433, codes['Sa_NCh433'][1,2], label='NCh433+DS61 - Zona 2')

axb2.semilogx(Tn_2745, codes['Sa_NCh2745'][1,0], label='NCh2745 - Zona 2')
axc2.semilogx(Tn_2745, codes['Sa_NCh2745'][1,1], label='NCh2745 - Zona 2')
axd2.semilogx(Tn_2745, codes['Sa_NCh2745'][1,2], label='NCh2745 - Zona 2')

axb2.text(0.05, 0.8, 'Suelo Tipo B', transform = axb2.transAxes)
axc2.text(0.05, 0.8, 'Suelo Tipo C', transform = axc2.transAxes)
axd2.text(0.05, 0.8, 'Suelo Tipo D', transform = axd2.transAxes)

axb2.set_ylabel('Sa [g]')
axc2.set_ylabel('Sa [g]')
axd2.set_ylabel('Sa [g]')

axb2.legend(loc=1)
axc2.legend(loc=1)
axd2.legend(loc=1)

axb2.set_xlim([1e-1, 10])
axc2.set_xlim([1e-1, 10])
axd2.set_xlim([1e-1, 10])

#axb2.set_ylim([0, 5.13])
#axc2.set_ylim([0, 5.13])
#axd2.set_ylim([0, 5.13])

axd2.set_xlabel('Periodo [s]')
fig2.suptitle('Aceleración espectral RotD100')

axb3.semilogx(Tn_433[:pos], codes['Sd_NCh433'][1,0][:pos]/100, label='NCh433+DS61 - Zona 2')
axc3.semilogx(Tn_433[:pos], codes['Sd_NCh433'][1,1][:pos]/100, label='NCh433+DS61 - Zona 2')
axd3.semilogx(Tn_433[:pos], codes['Sd_NCh433'][1,2][:pos]/100, label='NCh433+DS61 - Zona 2')

axb3.semilogx(Tn_433[pos:], codes['Sd_NCh433'][1,0][pos:]/100, label='', c=list(axb3.get_lines())[-1].get_color(), ls='--')
axc3.semilogx(Tn_433[pos:], codes['Sd_NCh433'][1,1][pos:]/100, label='', c=list(axc3.get_lines())[-1].get_color(), ls='--')
axd3.semilogx(Tn_433[pos:], codes['Sd_NCh433'][1,2][pos:]/100, label='', c=list(axd3.get_lines())[-1].get_color(), ls='--')

axb3.semilogx(Tn_2745, codes['Sd_NCh2745'][1,0]/100, label='NCh2745 - Zona 2')
axc3.semilogx(Tn_2745, codes['Sd_NCh2745'][1,1]/100, label='NCh2745 - Zona 2')
axd3.semilogx(Tn_2745, codes['Sd_NCh2745'][1,2]/100, label='NCh2745 - Zona 2')

axb3.text(0.8, 0.8, 'Suelo Tipo B', transform = axb3.transAxes)
axc3.text(0.8, 0.8, 'Suelo Tipo C', transform = axc3.transAxes)
axd3.text(0.8, 0.8, 'Suelo Tipo D', transform = axd3.transAxes)

axb3.set_ylabel('Sd [m]')
axc3.set_ylabel('Sd [m]')
axd3.set_ylabel('Sd [m]')

axb3.legend(loc=2)
axc3.legend(loc=2)
axd3.legend(loc=2)

axb3.set_xlim([1e-1, 10])
axc3.set_xlim([1e-1, 10])
axd3.set_xlim([1e-1, 10])

#axb3.set_ylim([0, 4.3])
#axc3.set_ylim([0, 4.3])
#axd3.set_ylim([0, 4.3])

axd3.set_xlabel('Periodo [s]')
fig3.suptitle('Desplazamiento espectral RotD50')

axb4.semilogx(Tn_433[:pos], codes['Sd_NCh433'][1,0][:pos]/100, label='NCh433+DS61 - Zona 2')
axc4.semilogx(Tn_433[:pos], codes['Sd_NCh433'][1,1][:pos]/100, label='NCh433+DS61 - Zona 2')
axd4.semilogx(Tn_433[:pos], codes['Sd_NCh433'][1,2][:pos]/100, label='NCh433+DS61 - Zona 2')

axb4.semilogx(Tn_433[pos:], codes['Sd_NCh433'][1,0][pos:]/100, label='', c=list(axb4.get_lines())[-1].get_color(), ls='--')
axc4.semilogx(Tn_433[pos:], codes['Sd_NCh433'][1,1][pos:]/100, label='', c=list(axc4.get_lines())[-1].get_color(), ls='--')
axd4.semilogx(Tn_433[pos:], codes['Sd_NCh433'][1,2][pos:]/100, label='', c=list(axd4.get_lines())[-1].get_color(), ls='--')

axb4.semilogx(Tn_2745, codes['Sd_NCh2745'][1,0]/100, label='NCh2745 - Zona 2')
axc4.semilogx(Tn_2745, codes['Sd_NCh2745'][1,1]/100, label='NCh2745 - Zona 2')
axd4.semilogx(Tn_2745, codes['Sd_NCh2745'][1,2]/100, label='NCh2745 - Zona 2')

axb4.text(0.8, 0.8, 'Suelo Tipo B', transform = axb4.transAxes)
axc4.text(0.8, 0.8, 'Suelo Tipo C', transform = axc4.transAxes)
axd4.text(0.8, 0.8, 'Suelo Tipo D', transform = axd4.transAxes)

axb4.set_ylabel('Sd [m]')
axc4.set_ylabel('Sd [m]')
axd4.set_ylabel('Sd [m]')

axb4.legend(loc=2)
axc4.legend(loc=2)
axd4.legend(loc=2)

axb4.set_xlim([1e-1, 10])
axc4.set_xlim([1e-1, 10])
axd4.set_xlim([1e-1, 10])

#axb4.set_ylim([0, 4.3])
#axc4.set_ylim([0, 4.3])
#axd4.set_ylim([0, 4.3])

axd4.set_xlabel('Periodo [s]')
fig4.suptitle('Desplazamiento espectral RotD100')

axb5.semilogx(Tn, ratio_b.mean(0), c='k', label='Media')
axb5.semilogx([Tn[0], Tn[-1]], [1.2, 1.2], c='tab:blue', ls='--', label='')
axb5.semilogx([Tn[0], Tn[-1]], [1.3, 1.3], c='tab:blue', ls='--', label='')

axc5.semilogx(Tn, ratio_c.mean(0), c='k', label='Media')
axc5.semilogx([Tn[0], Tn[-1]], [1.2, 1.2], c='tab:blue', ls='--', label='')
axc5.semilogx([Tn[0], Tn[-1]], [1.3, 1.3], c='tab:blue', ls='--', label='')

axd5.semilogx(Tn, ratio_d.mean(0), c='k', label='Media')
axd5.semilogx([Tn[0], Tn[-1]], [1.2, 1.2], c='tab:blue', ls='--', label='')
axd5.semilogx([Tn[0], Tn[-1]], [1.3, 1.3], c='tab:blue', ls='--', label='')


axb5.text(0.05, 0.8, 'Suelo Tipo B', transform = axb5.transAxes)
axc5.text(0.05, 0.8, 'Suelo Tipo C', transform = axc5.transAxes)
axd5.text(0.05, 0.8, 'Suelo Tipo D', transform = axd5.transAxes)

axb5.set_ylabel('RotD100/RotD50')
axc5.set_ylabel('RotD100/RotD50')
axd5.set_ylabel('RotD100/RotD50')

axb5.legend(loc=1)
axc5.legend(loc=1)
axd5.legend(loc=1)

axb5.set_xlim([1e-1, 10])
axc5.set_xlim([1e-1, 10])
axd5.set_xlim([1e-1, 10])

#axb1.set_ylim([0, 5.13])
#axc1.set_ylim([0, 5.13])
#axd1.set_ylim([0, 5.13])

axd5.set_xlabel('Periodo [s]')
fig5.suptitle('Cuociente RotD100/RotD50')


fig1.savefig('rotd50_2010_sa.pdf', bbox_inches='tight', pad_inches=0)
fig2.savefig('rotd100_2010_sa.pdf', bbox_inches='tight', pad_inches=0)
fig3.savefig('rotd50_2010_sd.pdf', bbox_inches='tight', pad_inches=0)
fig4.savefig('rotd100_2010_sd.pdf', bbox_inches='tight', pad_inches=0)
fig5.savefig('cuo_rotd100_rotd50_2010_sd.pdf', bbox_inches='tight', pad_inches=0)
