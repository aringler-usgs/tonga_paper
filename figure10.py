#!/usr/bin/env python
from obspy.core import UTCDateTime, read, Stream
import matplotlib.pyplot as plt
from scipy.signal import periodogram
import numpy as np
import math
from obspy.clients.fdsn import Client
from scipy import signal

import matplotlib as mpl
mpl.rc('font', family='serif')
mpl.rc('font', serif='Times')
mpl.rc('text', usetex=True)
mpl.rc('font', size=18)



mf, Mf = 2, 6
hours = 26
stime = UTCDateTime(UTCDateTime('2022-01-15T04:14:00'))


st = read('./tonga_bp.mseed')


st.detrend('linear')
st.merge(fill_value=0)
st.trim(stime,stime+hours*60*60)
for tr in st:
    tr.data *= 100


fig, ax = plt.subplots(2, 1, figsize=(16, 12))
ax = ax.flatten()
st2 = st.copy()


# filter between 2.5 and 5 mHz
st2.filter('bandpass', freqmin=1/400, freqmax=1/200)
ax[0].plot(st2[0].times()/(60*60), st2[0].data,
                 alpha=0.5, color='C0', label='AGBOM TONGA')


f_mHz = 3.68

days = 1.5
wo = 2*np.pi*(f_mHz/1000.)
Ao = np.max(st2[0].data)/2.
Q_True = 117

#t = np.arange(0,24*days*3600.)
t = st2[0].times()

Synthetic = Ao*np.exp((-wo*t)/(2*Q_True))*np.sin(wo*t)


ax[0].plot(st2[0].times()/(60*60), Synthetic,
                 alpha=0.5, color='C8', label='Synthetic Q=117')
ax[0].set_ylabel('Pressure (Pa)')
ax[0].set_xlabel('Time (hr.)')
ax[0].set_xlim((min(st2[0].times()/(60*60)), max(st2[0].times()/(60*60))))
ax[0].legend(loc='upper right', ncol=2)


# Window
for tr in st:
    tr.data *= signal.get_window(('kaiser', 2. * np.pi), tr.stats.npts)

NFFT = 2 ** (math.ceil(math.log(st[0].stats.npts, 2)))
print('Window Length is:', NFFT)


for idx, tr in enumerate(st):

    f, p = periodogram(tr.data, fs=tr.stats.sampling_rate,
                       nfft=NFFT, scaling='spectrum')
    p, f = p[1:], f[1:]

    # Now have p in nm/s/s switch f to mHz
    f *= 1000.
    p = p[(f >= mf) & (f <= Mf)]
    f = f[(f >= mf) & (f <= Mf)]
    if idx == 0:
        ln1 = ax[1].plot(f, p, label='AGBOM TONGA',
                         alpha=0.3, color='C0')
        ax[1].fill_between(f, 0, p, alpha=0.7)
        ax[1].set_xlim((mf, Mf))
        ax[1].set_ylabel('Amplitude (Pa)')
        #ax[1].tick_params(axis='y')
        ax[1].set_xlabel('Frequency (mHz)')
        ax[1].set_ylim((0, 1.1*max(p)))
    else:
        ax2 = ax[1].twinx()
        ln2 = ax2.plot(f, p, label='AGBOM TONGA',
                       alpha=0.3, color='C5')
        ax2.fill_between(f, 0, p, alpha=0.7, color='C5')
        ax2.set_ylabel('Amplitude ($\mu m/s)$', color='C5')
        ax2.tick_params(axis='y', labelcolor='C5')
        ax2.set_ylim((0, 1.1*max(p)))

# Modes from PREM calculated by MINEOS
ax[1].plot([3.68, 3.68], [-1, 40], color='C1', alpha=0.7)
ax[1].plot([4.40, 4.40], [-1, 40], color='C2', alpha=0.7)
ax[1].plot([3.63, 3.63], [-1, 40], color='C3', alpha=0.7)
ax[1].plot([3.72, 3.72], [-1, 40], color='C4', alpha=0.7)


ax[1].text(3.68, 0.3, '$_0S_{28} - _0S_{29}$',
           color='C1', alpha=0.7, ha='center', fontsize=22)
ax[1].text(4.40, 0.3, '$_0S_{36} - _0S_{37}$',
           color='C2', alpha=0.7, ha='center', fontsize=22)
ax[1].text(3.53, .25, '$_0S_{28}$', color='C3',
           alpha=0.7, ha='center', fontsize=22)
ax[1].text(3.82, .25, '$_0S_{29}$', color='C4',
           alpha=0.7, ha='center', fontsize=22)
ax[1].text(1.62, 0.3, '(b)')
ax[0].text(-2.5, 410, '(a)')

plt.savefig('Figure10.PNG', format='PNG', dpi=400)
plt.savefig('Figure10.PDF', format='PDF', dpi=400)

#plt.show()
