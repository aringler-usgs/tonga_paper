#!/usr/bin/env python
from obspy.core import read, UTCDateTime, Stream
import matplotlib.pyplot as plt
from scipy.signal import periodogram
from obspy.signal.invsim import evalresp
import numpy as np
import math
from obspy.clients.fdsn import Client 
from scipy import signal
import cartopy.crs as ccrs
import cartopy as cart
debug = True

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)


# network
client = Client("IRIS")
eve = UTCDateTime('2022-01-15T04:14:00')
chan_code = "LHZ"



mf, Mf= 3.5, 3.8
hours = 12


inv = client.get_stations(network='IU,IC,CU,II', station='*', starttime=eve, 
                endtime = eve + hours*60*60, level="response", location='*', channel='LHZ,LDI,LDO')

gidx = 0
idxall =0
lats, lons, press, hasdata = [], [], [], []
for net in inv:
    for sta in net:
        idxall+=1
        lats.append(sta.latitude)
        lons.append(sta.longitude)
        haspress = False
        for chan in sta:
            if chan.code in ['LDO','LDI']:
                haspress = True
        press.append(haspress)
        try:
            st = client.get_waveforms(network='IU,IC,CU,II', station=sta.code, location='00',
                                  channel='LHZ', starttime=eve,
                                  endtime=eve + hours*60*60)
            hasdata.append(True)
            gidx +=1
        except:
            hasdata.append(False)

print(str(gidx) + ' ' + str(idxall))
fig = plt.figure(1, figsize=(12,8))
ax = fig.add_subplot(1,1,1, projection = ccrs.Robinson())
ax.coastlines()
ax.add_feature(cart.feature.LAND, zorder=2, edgecolor='k')
ax.set_global()
ax.coastlines()
ax.scatter(175.4, -20.5, c='r',marker='*',s= 200, transform=ccrs.Geodetic(), zorder=3, vmin=0., vmax=10, label='Hunga Tonga Eruption')
ax.scatter(180-175.4, 20.5, c='k',marker='*',s= 200, transform=ccrs.Geodetic(), zorder=3, vmin=0., vmax=10, label='Eruption Antipode')
beenherepress = False
beenhereseis = False
beenheredata = False
presscount = 0
for trip in zip(lats, lons, press, hasdata):
    lat, lon, cpress, data = trip[0], trip[1], trip[2], trip[3]
    if data:

        if cpress:
            presscount +=1
            if beenherepress:
                ax.scatter(lon, lat, c='C0', s= 200, transform=ccrs.Geodetic(), zorder=3, alpha=0.5)

            else:
                ax.scatter(lon, lat, c='C0', s= 200, transform=ccrs.Geodetic(), zorder=3, alpha=0.5, label='Pressure and Seismic')
                beenherepress = True
        else:
            if beenhereseis:
                ax.scatter(lon, lat, c='C1', s= 200, transform=ccrs.Geodetic(), zorder=3, alpha=0.5)
            else:
                ax.scatter(lon, lat, c='C1', s= 200, transform=ccrs.Geodetic(), zorder=3, alpha=0.5, label='Seismic')
                beenhereseis = True
    else:
        if beenheredata:
            ax.scatter(lon, lat, c='C2', s = 200, transform=ccrs.Geodetic(), zorder=3, alpha=0.5)

        else:
            ax.scatter(lon, lat, c='C2', s = 200, transform=ccrs.Geodetic(), zorder=3, alpha=0.5, label='No Data')
            beenheredata = True
fig.legend(ncol=5, fontsize=14, loc='lower center')
plt.savefig('Figure1.PNG', format='PNG', dpi=400)
plt.savefig('Figure1.PDF', format='PDF', dpi=400)

print('Here are the pressure stations:' + str(presscount))
