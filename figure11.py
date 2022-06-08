#!/usr/bin/env python
from obspy import UTCDateTime, read_inventory, read
import matplotlib.pyplot as plt
from scipy.signal import periodogram
import numpy as np
import math
from obspy.clients.fdsn import Client
from scipy import signal

import matplotlib as mpl
mpl.rc('font', family='serif')
mpl.rc('font', serif='Times')
#mpl.rc('text', usetex=True)
mpl.rc('font', size=18)

debug = True

'''
Request 00 location code LDO and LHZ data from II MSVF for 26 hours
starting at the event start time.


remove the linear trend
zero fill data gaps
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


'''
Minimum and Maximum frequencies to plot
Hours to calculate spectra for following Event_time
'''

mf, Mf = 3.5, 4
hours = 26
FR_trim = 0

# Manually pick start index for this one off (row of CSV corresponding to stime)
SI = 254

# Get AFI Data

inv = client.get_stations(network='IU,IC,CU,II', station='MSVF', starttime=eve,
                          endtime=eve + hours*60*60, level="response",
                          location='00', channel='LDO')

st = client.get_waveforms(network='IU,IC,CU,II', station='MSVF', location='00',
                          channel='LDO', starttime=eve - FR_trim*60*60,
                          endtime=eve + hours*60*60)


st.detrend('linear')
st.merge(fill_value=0)


# Now we get the tonga data

st_Tonga = read('./tonga_bp.mseed')

st_Tonga.detrend('linear')
st_Tonga.merge(fill_value=0)
st_Tonga.trim(eve-FR_trim*60*60,eve+hours*60*60)

if debug:
    st_Tonga.plot()

# Convert from mbar to Pa
st_Tonga[0].data *= 100


fig, ax = plt.subplots(1, 1, figsize=(16, 12))


# Window
for tr in st:
    tr.data *= signal.get_window(('kaiser', 2. * np.pi), tr.stats.npts)

for tr in st_Tonga:
    tr.data *= signal.get_window(('kaiser', 2. * np.pi), tr.stats.npts)


NFFT = 2 ** (math.ceil(math.log(st[0].stats.npts, 2)))
NFFT_Tonga = 2 ** (math.ceil(math.log(st_Tonga[0].stats.npts, 2)))

print('Window Length is:', NFFT)

st.append(st_Tonga[0])

print(st)

colors = ['C0', 'C1']

for idx, tr in enumerate(st):

    f, p = periodogram(tr.data, fs=tr.stats.sampling_rate,
                       nfft=NFFT, scaling='spectrum')
    p, f = p[1:], f[1:]

    # If MSVF we, remove response by spectral division
    if tr.stats.station == 'MSVF':
        print('Found MSVF')
        trid = tr.id
        inv_resp = inv.get_response(trid, tr.stats.starttime)
        resp, _ = inv_resp.get_evalresp_response(tr.stats.delta, NFFT)
        resp = resp[1:]
        p = np.sqrt(p/(np.abs(resp)**2))
    # If Tonga data, we already have the data in Pa
    else:
        p = np.sqrt(p)

    # Now have p in nm/s/s switch f to mHz
    f *= 1000.
    p = p[(f >= mf) & (f <= Mf)]
    f = f[(f >= mf) & (f <= Mf)]
    if idx == 0:

        ax.plot(f, p, label=(tr.id).replace('.', ' '),
                         alpha=0.3, color=colors[idx])
    else:
        ax.plot(f, p, label='Tonga Barometer',
                         alpha=0.3, color=colors[idx])

    ax.fill_between(f, 0, p, alpha=0.7, color=colors[idx])
    ax.set_xlim((mf, Mf))
    ax.set_ylabel('Amplitude (Pa)', color='C0')
    ax.tick_params(axis='y', labelcolor='C0')
    ax.set_xlabel('Frequency (mHz)')
    ax.set_ylim((0, 1.1*max(p)))


# Modes from PREM calculated by MINEOS
ax.plot([3.68, 3.68], [-1, 2], color='C1', alpha=0.7)
ax.plot([3.63, 3.63], [-1, 2], color='C3', alpha=0.7)
ax.plot([3.72, 3.72], [-1, 2], color='C4', alpha=0.7)


ax.legend(loc=0)
ax.text(3.68, 0.55, '$_0S_{28} - _0S_{29}$',
           color='C1', alpha=0.7, ha='center', fontsize=22)
ax.text(3.614, 0.54, '$_0S_{28}$', color='C3',
           alpha=0.7, ha='center', fontsize=22)
ax.text(3.705, 0.585, '$_0S_{29}$', color='C4',
           alpha=0.7, ha='center', fontsize=22)

plt.savefig('Figure11.PNG', format='PNG', dpi=400)
plt.savefig('Figure11.PDF', format='PDF', dpi=400)

#plt.show()
