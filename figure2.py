#!/usr/bin/env python
from obspy.core import UTCDateTime, AttribDict
from obspy.clients.fdsn import Client 
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

client = Client("IRIS")
nets = 'IU,II,IC,CU'
stime = UTCDateTime('2022-01-15T00:00:00')
ev_coord = [-20.5, 175.4, ]

chans, locs = ['LDO,LDI','LHZ'], ['31,00', '00']

fig, ax = plt.subplots(1, 2, figsize=(16,16))
ax = ax.flatten()
for idx, chanloc in enumerate(zip(chans, locs)):
    chan, loc = chanloc
    inv = client.get_stations(network=nets, starttime=stime, 
            endtime = stime + 60*60*60, level="response", 
            location=loc, channel=chan, station='*')

    st = client.get_waveforms(network=nets, starttime=stime, 
        endtime = stime + 60*60*60, 
        location=loc, channel=chan,station='*')
    st.detrend('constant')
    st.merge(fill_value=0.)
    if chan == 'LHZ':
        st.remove_response(inventory=inv)
    st.taper(0.05)
    st.filter('bandpass',freqmin=0.001, freqmax=0.1)
    for tr in st:
        coors = inv.get_coordinates(tr.id, tr.stats.starttime)
        tr.stats.coordinates = AttribDict()
        tr.stats.coordinates.latitude = coors['latitude']
        tr.stats.coordinates.longitude = coors['longitude']

    plt.axes(ax[idx])
    st.plot(type='section', ev_coord=ev_coord, show=False, dist_degree=True,orientation='horizontal',fig=fig,alph=0.5)
    ax[idx].set_title('')

    if idx == 0:
        ax[idx].set_ylabel('Offset ($\deg$)', fontsize=18)
    else:
        ax[idx].set_ylabel('')
    
    ax[idx].set_xlabel('Time (hr.)', fontsize=18)
    #ax[idx].set_xticks([1,2])
    vals = ax[idx].get_xticklabels()
    for val in vals:
        cval = val.get_text()
        cval = cval.replace('$\mathdefault{','')
        cval = cval.replace('}$','')
        cval = str(np.round(float(cval)/(60*60),2))
        val.set_text('$\mathdefault{' + cval + '}$')
    ax[idx].set_xticklabels(vals)

ax[0].text(-14100, 170,'(a)')
ax[1].text(-14100, 170, '(b)')    







plt.savefig('Figure2.png',format='PNG', dpi=400)
plt.savefig('Figure2.pdf',format='PDF', dpi=400)