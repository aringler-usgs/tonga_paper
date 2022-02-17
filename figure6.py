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
eve = UTCDateTime('2022-01-15T04:14:00')
chan_code = "LDO,LDI"



mf, Mf= 3, 6
hours = 26


inv = client.get_stations(network='IU,IC,CU,II', station='*', starttime=eve, 
                endtime = eve + hours*60*60, level="response", location='00,31', channel=chan_code)


st = client.get_waveforms(network='IU,IC,CU,II', station='*', location='00,31', channel=chan_code, 
            starttime=eve, endtime=eve + hours*60*60)

st.merge(fill_value=0)
st.detrend('linear')
st.detrend('constant')
st.sort(['station'])

fig, ax = plt.subplots(len(st), 1, figsize=(16,12))
ax = ax.flatten()

# Window
for tr in st:
    tr.data *= signal.get_window(('kaiser', 2.*np.pi), tr.stats.npts)

NFFT=2**(math.ceil(math.log(st[0].stats.npts, 2)))

fhand = open('PressureSNR.csv','w')

for idx, tr in enumerate(st):
    

    f,p = periodogram(tr.data, fs=tr.stats.sampling_rate, nfft= NFFT, scaling='spectrum')
    p, f = p[1:], f[1:]
    trid = tr.id
    #inv_resp = inv.get_response(trid, tr.stats.starttime)
    #resp, _ = inv_resp.get_evalresp_response(tr.stats.delta, NFFT, 'ACC')
    #resp = resp[1:]
    # Convert units to nm/s/s
    p = np.sqrt(p)*10**9
    #p = np.sqrt(p/(np.abs(resp)**2))*10**9
    # Now have p in nm/s/s switch f to mHz
    f *= 1000.
    p = p[(f >= mf) & (f <= Mf)]
    f = f[(f >= mf) & (f <= Mf)]
    snr1 = np.max(p[(f >= 3.5) & (f <= 3.8)])/np.max(p[(f >= 4.0) & (f <= 4.3)])
    snr2 = np.max(p[(f >= 4.3) & (f <= 4.6)])/np.max(p[(f >= 5.0) & (f <= 5.3)])

    if 'MSVF' in trid:
        print(trid)
        print(snr1)
        print(snr2)
    if snr1 > 3:
        print(snr1)
        print(trid)



    fhand.write(trid + ', ' + str(snr1) + ', ' + str(snr2) + '\n')


    p /= np.max(p)
    ax[idx].plot(f,p, label=(tr.id).replace('.',' '), alpha=0.3, color='C0')
    ax[idx].fill_between(f, 0, p, alpha=0.7, color='C0')
    ax[idx].set_xlim((mf,Mf))
    ax[idx].plot([3.68, 3.68], [-1,2],color='C1',alpha=0.7)
    ax[idx].plot([4.40, 4.40], [-1,2],color='C2',alpha=0.7)
    ax[idx].plot([3.63, 3.63], [-1,2],color='C3',alpha=0.7)
    ax[idx].plot([3.72, 3.72], [-1,2],color='C4',alpha=0.7)
    ax[idx].set_ylim((0,1.1))
    ax[idx].set_yticks([0.55])
    ax[idx].set_yticklabels([tr.stats.station], fontsize=8)

    if idx < len(st)-1:
        ax[idx].set_xticks([])
#ax.set_yticks([0.5])
fhand.close()

ax[0].text(3.68, 9, '$_0S_{28} - _0S_{29}$', color='C1', alpha=0.7, ha='center')
ax[0].text(4.40, 3, '$_0S_{36} - _0S_{37}$', color='C2', alpha=0.7, ha='center')
ax[0].text(3.63, 5, '$_0S_{28}$', color='C3', alpha=0.7, ha='center')
ax[0].text(3.72, 3, '$_0S_{29}$', color='C4', alpha=0.7, ha='center')
ax[idx].set_xlabel('Frequency ($mHz$)')
plt.savefig('Figure6.PNG', format='PNG', dpi=400)
plt.savefig('Figure6.PDF', format='PDF', dpi=400)
plt.close('all') 

    #plt.savefig('Test' + '.pdf', format='PDF', dpi=400) 