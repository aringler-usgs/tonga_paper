#!/usr/bin/env python
from obspy.core import UTCDateTime, AttribDict
from obspy.clients.fdsn import Client 
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from obspy.geodetics import kilometers2degrees, gps2dist_azimuth

mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

client = Client("IRIS")
nets = 'CU,II,IC,IU'
stime = UTCDateTime('2022-01-15T00:00:00')
ev_coord = [-20.5, 175.4, ]

chans, locs = ['LDO,LDI','LHZ'], ['31,00', '00']

lamb = kilometers2degrees(0.314)*60*60
surface = kilometers2degrees(4)*60*60

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
    if idx == 0:
        st.filter('bandpass',freqmin=0.001, freqmax=0.1, zerophase=True)
    else:
        st.filter('bandpass',freqmin=0.001, freqmax=0.1, zerophase=True)
    for tr in st:
        coors = inv.get_coordinates(tr.id, tr.stats.starttime)
        tr.stats.coordinates = AttribDict()
        tr.stats.coordinates.latitude = coors['latitude']
        tr.stats.coordinates.longitude = coors['longitude']
        dist,_, _ = gps2dist_azimuth(coors['latitude'], coors['longitude'], ev_coord[0], ev_coord[1])
        dist = kilometers2degrees(dist/1000.)
        #ax[idx].plot(tr.times()/(60*60), 3*tr.data/np.abs(np.max(tr.data)) + dist, color='0.4', alpha=0.7)
    plt.axes(ax[idx])
    



    axes = st.plot(type='section', ev_coord=ev_coord, show=False, dist_degree=True,orientation='horizontal',fig=fig,alph=0.5)
    plt.plot([4.15*60*60, (17.93 + 4.15)*60*60], [0, 17.93*lamb],color='C0')
    #plt.plot([(17.93 + 4.15)*60*60, (2*18.93 + 4.15)*60*60], [180., 0.],color='C0')
    #ax[idx].plot([18+4.15, 18+4.15], [0,180],color='C2')
    plt.plot([4.15*60*60,(20+ 4.15)*60*60], [0, (20 + 4.15)*surface],color='C1')
    #print(axes)
    ax[idx].set_title('')

    if idx == 0:
        ax[idx].set_ylabel('Offset ($\deg$)', fontsize=18)
    else:
        ax[idx].set_ylabel('')
    
    ax[idx].set_xlabel('Time (hr.)', fontsize=18)
    #ax[idx].set_xticks([1,2])
    vals = ax[idx].get_xticklabels()
    from matplotlib.text import Text 
    vals = [Text(0.0, 0., '$\\mathdefault{0}$'), 
        Text(10.0, 0., '$\\mathdefault{10}$'),Text(20., 0., '$\\mathdefault{20}$'),
        Text(30.0, 0., '$\\mathdefault{30}$'),Text(40.0, 0., '$\\mathdefault{40}$'),Text(50, 0., '$\\mathdefault{50}$')]


    #for val in vals:
    #    cval = val.get_text()
    #    cval = cval.replace('$\mathdefault{','')
    #    cval = cval.replace('}$','')
    #    cval = str(np.round(float(cval)/(60*60),2))
    #    val.set_text('$\mathdefault{' + cval + '}$')
    ax[idx].set_xticklabels(vals)
    ax[idx].set_ylim((0,180))
    #ax[idx].set_xlim((0,50))
ax[0].text(-14100, 183,'(a)')
ax[1].text(-14100, 183, '(b)')    







plt.savefig('Figure2.png',format='PNG', dpi=400)
plt.savefig('Figure2.pdf',format='PDF', dpi=400)

#plt.show()