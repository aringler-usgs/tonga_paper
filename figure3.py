#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import obspy

from obspy import UTCDateTime
from obspy.clients.fdsn import Client 

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

#Global Variables
# Distance from ANMO to Tonga (km)

debug = True

ANMO_dist = 9524
eve = UTCDateTime('2022-015T04:14:45')
stime = UTCDateTime('2022-015T00:00:00')
# Rayleigh wave velocity and Acoustic velocity (km/s)
# Pressure taken from Junghyun Park, SMU estimate

R_speed = 4.0
A_speed = 0.32
# P-wave arrival time delay (s)
P_delay = 762

# Plot delay to remove filter ringing - 2 hours
Plot_delay = 2.*3600.

##############################################################################


client = Client("IRIS")

inv = client.get_stations(network='IU', station='ANMO', starttime=stime, 
                endtime = stime + 24*60*60, level="response", location='*', channel='LHZ,LDI,LDO')




st = client.get_waveforms(network='IU', station='ANMO', location='*',
                                  channel='LHZ,LDI,LDO', starttime=stime,
                                  endtime=stime + 24*60*60)



st.detrend('constant')
st.remove_response(inventory=inv)

st2 = st.copy()
st3 = st.copy()

st.filter('bandpass',freqmin=0.003, freqmax = 0.05 ,corners=4)
st2.filter('bandpass',freqmin=0.01, freqmax = 0.04 ,corners=4)
st3.filter('bandpass',freqmin=0.003, freqmax = 0.005 ,corners=4)

st.trim(st[0].stats.starttime + Plot_delay,st[0].stats.endtime)
st2.trim(st[0].stats.starttime,st[0].stats.endtime)
st3.trim(st[0].stats.starttime, st[0].stats.endtime)

# Convert seismic data to micrometers per second
st.select(channel='LHZ')[0].data *= 1000000
st2.select(channel='LHZ')[0].data *= 1000000
st3.select(channel='LHZ')[0].data *= 1000000


######## Make the Figure ###########################################
fig = plt.figure(1, figsize=(12,16))
#plt.subplots_adjust(hspace=0.001)

color1 = 'C0'
ax1 = plt.subplot(211)
ax1.plot((st.select(channel = 'LDO')[0].times()/3600.)+2., st.select(channel = 'LDO')[0].data,linewidth=1.5,c=color1,alpha=1.0)
ax1.set_ylabel('Pressure (Pa)', fontsize=18, c=color1)
ax1.tick_params(axis='y', colors=color1, labelsize=15)
ax1.set_ylim(-30., 30.)
ax1.set_xlim(2.0,24.0)
ax1.set_xlabel('Time (hr.)', fontsize=18)
ax1.text(0.0, 30, '(a)')
# since we have a bunch of plots with the same time, we only want
# xticks on the bottom
#plt.tick_params(labelbottom=False)


ax2 = plt.subplot(212)
ax2.plot((st.select(channel = 'LHZ')[0].times()/3600.)+2, st.select(channel = 'LHZ')[0].data,linewidth=1.0,c='C1',alpha=1.0)
ax2.text(18,0.2, 'Bandpass 3 to 50 mHz', color='C1')
ax2.plot((st2.select(channel = 'LHZ')[0].times()/3600.)+2, st2.select(channel = 'LHZ')[0].data+1.,linewidth=1.0,c='C2',alpha=0.7)
ax2.text(18,1.2, 'Bandpass 10 to 40 mHz', color='C2')
ax2.plot((st3.select(channel = 'LHZ')[0].times()/3600.)+2, 2.0*(st3.select(channel = 'LHZ')[0].data)-1.,linewidth=1.0,c='C3',alpha=0.8)
ax2.text(18,-0.8, 'Bandpass 3 to 5 mHz', color='C3')
ax2.set_ylabel('Ground Velocity ($\mu$m/s)', fontsize=18)
plt.tick_params(labelsize = 15)
ax2.set_ylim(-1.5, 1.8)
ax2.set_xlim(2.0,24.0)
ax2.set_xlabel('Time (hr.)', fontsize=18)
ax2.text(0.0, 1.8, '(b)')

plt.savefig('Figure3.PDF', format='PDF', dpi=400)
plt.savefig('Figure3.PNG', format='PNG', dpi=400)

