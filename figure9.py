#!/usr/bin/env python
from obspy.core import UTCDateTime
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
mpl.rc('font', size=16)


'''
Request 00 location code LDO and LHZ data from II MSVF for 24 hours
starting at the event start time.

zero fill any gaps
remove the linear trend
remove the response to velocity and Pa
bandpass filter from 200 to 400 s period

For the spectra apply a 2pi Kaiser window
calculate DFT using the next power of 2 length NFFT
remove the response to units of nm/s^2 in acceleration
references frequencies are calculated from PREM using mineos
'''


# network
client = Client("IRIS")
eve = UTCDateTime('2022-01-15T04:14:00') 

mf, Mf = 2.0, 6
hours = 26

inv = client.get_stations(network='IU,IC,CU,II,G', station='MSVF', starttime=eve,
                          endtime=eve + hours*60*60, level="response",
                          location='00', channel='LHZ,LDI')



st = client.get_waveforms(network='IU,IC,CU,II,G', station='MSVF', location='00',
                          channel='LHZ,LDI', starttime=eve,
                          endtime=eve + hours*60*60)


st.merge(fill_value=0)
st.detrend('linear')
st.detrend('constant')
st.sort(['station'])


fig, ax = plt.subplots(3, 1, figsize=(16, 8))
ax = ax.flatten()
st2 = st.copy()

st2.remove_response(inventory=inv)
st2.filter('bandpass', freqmin=1/400, freqmax=1/200)

f_mHz = 3.68

days = 1.5
wo = 2*np.pi*(f_mHz/1000.)
Ao = np.max(st2[1].data)*10**6
Q_True = 186.9

#t = np.arange(0,24*days*3600.)
t = st2[1].times()

Synthetic = Ao*np.exp((-wo*t)/(2*Q_True))*np.sin(wo*t)





ax[1].plot(st2[0].times()/(60*60), st2[0].data,
                 alpha=0.7, color='C0', label=(st2[0].id).replace('.', ' '))

ax[2].plot(st2[1].times()/(60*60), Synthetic,
               alpha=0.5, color='C8', label='Synthetic Q=' + str(Q_True))
ax[2].plot(st2[1].times()/(60*60), st2[1].data*10**6,
               alpha=0.7, color='C5', label=(st2[1].id).replace('.', ' '))

ax[2].set_ylabel('Velocity ($\mu m/s$)')

ax[1].set_ylabel('Pressure (Pa)')

ax[1].set_xlabel('Time (hr.)')
ax[2].set_xlabel('Time (hr.)')
ax[1].set_xlim((min(st2[1].times()/(60*60)), max(st2[1].times()/(60*60))))
ax[2].set_xlim((min(st2[1].times()/(60*60)), max(st2[1].times()/(60*60))))
# ln = ln1 + ln2 + ln3
# labs = [ls.get_label() for ls in ln]
ax[1].legend(loc='upper right', ncol=2)
ax[2].legend(loc='upper right', ncol=2)

# Window
for tr in st:
    tr.data *= signal.get_window(('kaiser', 2. * np.pi), tr.stats.npts)

NFFT = 2 ** (math.ceil(math.log(st[0].stats.npts, 2)))


st[1].data *= 10**6
for idx, tr in enumerate(st):

    f, p = periodogram(tr.data, fs=tr.stats.sampling_rate,
                       nfft=NFFT, scaling='spectrum')
    p, f = p[1:], f[1:]

    inv_resp = inv.get_response(tr.id, tr.stats.starttime)
    resp, _ = inv_resp.get_evalresp_response(tr.stats.delta, NFFT)
    resp = resp[1:]
    p = np.sqrt(p/(np.abs(resp)**2))

    # Convert units to nm/s/s
    #p = np.sqrt(p)*10**9
    print(st)
    if tr.stats.channel == 'LDI':
        f2 = f*1000
        temp = p[(f2>=3.6) & (f2 <= 3.8)]
        temp2 =f2[(f2>=3.6) & (f2 <= 3.8)]
        for pair in zip(temp,temp2):
            print(pair)
    # Now have p in nm/s/s switch f to mHz
    f *= 1000.
    p = p[(f >= mf) & (f <= Mf)]
    f = f[(f >= mf) & (f <= Mf)]
    if idx == 0:
        ln1 = ax[0].plot(f, p, label=(tr.id).replace('.', ' '),
                         alpha=0.3, color='C0')
        ax[0].fill_between(f, 0, p, alpha=0.7, color='C0')
        ax[0].set_xlim((mf, Mf))
        ax[0].set_ylabel('Amplitude ($Pa^2$)', color='C0')
        ax[0].tick_params(axis='y', labelcolor='C0')
        ax[0].set_xlabel('Frequency (mHz)')
        ax[0].set_ylim((0, 1.1*max(p)))
    else:
        ax2 = ax[0].twinx()
        ln2 = ax2.plot(f, p, label=(tr.id).replace('.', ' '),
                       alpha=0.3, color='C5')
        ax2.fill_between(f, 0, p, alpha=0.7, color='C5')
        ax2.set_ylabel('Amplitude ($\mu m/s^2$)', color='C5')
        ax2.tick_params(axis='y', labelcolor='C5')
        ax2.set_ylim((0, 1.1*max(p)))

# Modes from PREM calculated by MINEOS
ax[0].plot([3.68, 3.68], [-1, 2], color='C1', alpha=0.7)
ax[0].plot([4.40, 4.40], [-1, 2], color='C2', alpha=0.7)
ax[0].plot([3.63, 3.63], [-1, 2], color='C3', alpha=0.7)
ax[0].plot([3.72, 3.72], [-1, 2], color='C4', alpha=0.7)

ln = ln1 + ln2
labs = [ls.get_label() for ls in ln]
ax2.legend(ln, labs, loc=0, ncol=2)
ax[0].text(3.68, 0.322+0.05, '$_0S_{28} - _0S_{29}$',
           color='C1', alpha=0.7, ha='center', fontsize=22)
ax[0].text(4.40, 0.322+0.05, '$_0S_{36} - _0S_{37}$',
           color='C2', alpha=0.7, ha='center', fontsize=22)
ax[0].text(3.63, 0.345+0.1, '$_0S_{28}$', color='C3',
           alpha=0.7, ha='center', fontsize=22)
ax[0].text(3.82, 0.342+0.1, '$_0S_{29}$', color='C4',
           alpha=0.7, ha='center', fontsize=22)
ax[0].text(1.62, 0.36, '(a)')
ax[1].text(-2.5, 140, '(b)')
ax[2].text(-2.5, 0.6, '(c)')


plt.tight_layout()
plt.savefig('Figure9.PNG', format='PNG', dpi=400)
plt.savefig('Figure9.PDF', format='PDF', dpi=400)
