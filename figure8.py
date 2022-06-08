#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import obspy
from obspy.io.xseed import Parser
from obspy import read_inventory, UTCDateTime, read, Stream


#Global Variables
# Distance from ANMO to Tonga (km)

debug = True
eve = UTCDateTime('2022-01-15T04:14:00')
Hours = 1.5;
##############################################################################

st = Stream()
# hacky since we only have the resps we want in this directory
inv = read_inventory('./RESP*')


# Load data
st += read("*.mseed")
st.detrend('constant')

# Remove Response to Displacement and Pressure
st.select(channel='LHZ')[0].remove_response(inventory=inv,output='ACC')
st.select(channel='LDI')[0].remove_response(inventory=inv)

st.filter('bandpass',freqmin=0.003, freqmax = 0.005 ,corners=4)

st.trim(eve,eve+Hours*60*60)
# Convert seismic data (displacement) to nm.s/s
st.select(channel='LHZ')[0].data *= 10**9
st.select(channel='LDI')[0].data /= 100.

def newfun(a):
    return np.sum((st.select(channel='LHZ')[0].data - a*st.select(channel='LDI')[0].data)**2)

from scipy.optimize import fmin 

out = fmin(newfun, 3.5)
print(out)

print('Here we go')
print(newfun(out[0]))
print(newfun(0))
print(newfun(out[0])/newfun(0))


######## Make the Figure ###########################################
fig = plt.figure(1, figsize=(12,16))
ax1 = plt.subplot(111)
ax1.plot(st.select(channel='LHZ')[0].times()/60.,st.select(channel='LHZ')[0].data,linewidth=3.,label = 'II MSVF 00 LHZ', color='C1', alpha=0.7)
ax1.plot(st.select(channel='LHZ')[0].times()/60.,st.select(channel='LHZ')[0].data - out[0]*st.select(channel='LDI')[0].data,linewidth=1.,label = 'II MSVF 00 LHZ Pressure Corrected', color='C2', alpha=0.8)

# Quick hack to deal with twin X and legend
ax1.plot([], [], linewidth=1.5, color='C0', label = 'II MSVF 00 LDI')
ax1.set_ylabel('Ground\nAcceleration ($nm/s^2$)', fontsize=18)
plt.tick_params(labelsize = 15)
plt.xlabel('Time Since Eruption (min)', fontsize = 18)
plt.legend(fontsize=15)
ax2 = ax1.twinx()
ax2.plot(st.select(channel='LDI')[0].times()/60.,st.select(channel='LDI')[0].data,linewidth=3., color='C0',alpha=0.7)
ax2.set_ylabel('Pressure (hPa)', fontsize=18, c='C0')
ax2.tick_params(axis='y', color='C0', labelsize=15)
#ax2.set_ylim((-1.2,1.2))
plt.xlim((0, 1.5*60))
plt.legend()
plt.savefig('FigureS4_again.png', format='PNG', dpi=400)
plt.savefig('FigureS4_again.pdf', format='PDF', dpi=400)
#plt.show()
