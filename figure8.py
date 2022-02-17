#!/usr/bin/env python
from obspy.core import read, UTCDateTime, Stream
import matplotlib.pyplot as plt
from scipy.signal import periodogram, hilbert
from obspy.signal.invsim import evalresp
import numpy as np
from scipy.optimize import fmin
import math
import glob
import sys
from obspy.clients.fdsn import Client 
from scipy import signal
debug = True

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)




# network
client = Client("IRIS")
eve2 = UTCDateTime('2022-01-15T04:14:00')
eve1 = UTCDateTime('1991-06-15T06:30:00')
chan_code = "LHZ"



mf, Mf= 3, 5
hours = 12



inv = client.get_stations(network='II,IU,IC,CU', station='*', starttime=eve1, 
                endtime = eve1 + hours*60*60, level="response", location='*', channel=chan_code)

st1 = client.get_waveforms(network='II,IU,IC,CU', station='*', location='*', channel=chan_code, 
            starttime=eve1, endtime=eve1 + hours*60*60)

st1.merge(fill_value=0)
st1.detrend('linear')
st1.detrend('constant')


for tr in st1:
    if tr.stats.station == 'COL':
        st1.remove(tr)
    if tr.stats.station == 'GUMO':
        st1.remove(tr)

for tr in st1:
    tr.data *= signal.get_window(('kaiser', 2.*np.pi), tr.stats.npts)

NFFT=2**(math.ceil(math.log(st1[0].stats.npts, 2)))


for idx, tr in enumerate(st1):
    f,p = periodogram(st1[0].data, fs=tr.stats.sampling_rate, nfft= NFFT, scaling='spectrum')
    p, f = p[1:], f[1:]
    inv_resp = inv.get_response(tr.id, tr.stats.starttime)
    resp, _ = inv_resp.get_evalresp_response(tr.stats.delta, NFFT, 'ACC')
    resp = resp[1:]
    # Convert units to nm/s/s

    p = np.sqrt(p/(np.abs(resp)**2))*10**9
    # Now have p in nm/s/s switch f to mHz
    f *= 1000.
    p = p[(f >= mf) & (f <= Mf)]
    f = f[(f >= mf) & (f <= Mf)]
    #p /= np.max(p)
    #val = np.max([np.max(p), np.max(p2)])
    #rat = np.max(p)/np.max(p2)
    if idx == 0:
        avg1 = p
    else:
        avg1 += p 

avg1 /= (idx+1)


# Run code normally here to to plot each spectra.


inv = client.get_stations(network='II,IU,IC,CU', station='*', starttime=eve2, 
                endtime = eve2 + hours*60*60, level="response", location='*', channel=chan_code)


st = client.get_waveforms(network='II,IU,IC,CU', station='*', location='*', channel=chan_code, 
            starttime=eve2, endtime=eve2 + hours*60*60)



st.merge(fill_value=0)
st.detrend('linear')
st.detrend('constant')

fig, ax = plt.subplots(1, 1, figsize=(12,12))

for tr in st:
    tr.data *= signal.get_window(('kaiser', 2.*np.pi), tr.stats.npts)

NFFT=2**(math.ceil(math.log(st[0].stats.npts, 2)))


for idx, tr in enumerate(st):
    f,p = periodogram(tr.data, fs=tr.stats.sampling_rate, nfft= NFFT, scaling='spectrum')
    p, f = p[1:], f[1:]
    inv_resp = inv.get_response(tr.id, tr.stats.starttime)
    resp, _ = inv_resp.get_evalresp_response(tr.stats.delta, NFFT, 'ACC')
    resp = resp[1:]
    # Convert units to nm/s/s
    p = np.sqrt(p/(np.abs(resp)**2))*10**9
    # Now have p in nm/s/s switch f to mHz
    f *= 1000.
    p = p[(f >= mf) & (f <= Mf)]
    f = f[(f >= mf) & (f <= Mf)]
    #p /= np.max(p)

    check = p[(f>=3.9) & (f <=4.1)]
    check2 = p[(f>=3.0) & (f <=3.4)]

    if np.max(check) > 0.6 or np.max(check2) >0.35:
        continue
    if idx == 0:
        avg2 = p
    else:
        avg2 += p



    #p /= np.max(p)
    #p2 /= np.max(p2)



    if idx == 0:
        ax.plot(f,p,  color='C0', alpha=0.2, label='Spectra')
    else:
        ax.plot(f,p,  color='C0', alpha=0.2)
    ax.fill_between(f, 0, p, color='C0', alpha=0.2)
    
avg2 /= (idx+1)
ax.set_xlim((mf,Mf))
ax.plot([3.68, 3.68], [-1,2],color='C2',alpha=0.7)
ax.plot([4.40, 4.40], [-1,2],color='C3',alpha=0.7)
ax.plot([3.63, 3.63], [-1,2],color='C4',alpha=0.7)
ax.plot([3.72, 3.72], [-1,2],color='C5',alpha=0.7)
ax.plot([4.53, 4.53], [-1,2],color='C6',alpha=0.7)
ax.plot([4.62, 4.62], [-1,2],color='C7',alpha=0.7)
ax.set_ylim((0,1.5))


# p11 = np.max(avg1[(f >= 3.4) & (f <= 4.0)])
# p12 = np.max(avg1[(f >= 4.2) & (f <= 4.8)])
# p21 = np.max(avg2[(f >= 3.4) & (f <= 4.0)])
# p22 = np.max(avg2[(f >= 4.2) & (f <= 4.8)])


#ax.text(3.68+.1, p11+0.05, str(round(p11,2)) + ' $nm/s^2$', color='C1')
#ax.text(4.4+0.1, p12+0.05, str(round(p12,2))+ ' $nm/s^2$', color='C1')


ax.plot(f,avg1, color='C1', label='Pinatubo Mean', linewidth=3)
ax.plot(f,avg2, label='Tonga Mean', color='C8', linewidth=3)

ax.legend(loc ='upper left')
ax.text(3.68, 1.12+0.4, '$_0S_{28} - _0S_{29}$', color='C2', alpha=0.7, ha='center')
ax.text(4.40, 1.12+0.4, '$_0S_{36} - _0S_{37}$', color='C3', alpha=0.7, ha='center')
ax.text(3.63, 1.15+0.4, '$_0S_{28}$', color='C4', alpha=0.7, ha='center')
ax.text(3.72, 1.17+0.4, '$_0S_{29}$', color='C5', alpha=0.7, ha='center')

ax.text(4.53, 1.15+0.4, '$_0S_{38}$', color='C6', alpha=0.7, ha='center')
ax.text(4.62, 1.17+0.4, '$_0S_{39}$', color='C7', alpha=0.7, ha='center')
ax.set_xlabel('Frequency ($mHz$)')
ax.set_ylabel('Amplitude ($nm/s^2$)')

plt.savefig('Figure8.PNG', format='PNG', dpi=400)
plt.savefig('Figure8.PDF', format='PDF', dpi=400)


    #plt.savefig('Test' + '.pdf', format='PDF', dpi=400) 