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


def sph2cart(lat, lon, r):
    x = r * np.cos(lon*np.pi/180.) * np.cos(lat*np.pi/180.)
    y = r * np.sin(lon*np.pi/180.) * np.cos(lat*np.pi/180.)
    z = r * np.sin(lat*np.pi/180.)
    return [x, y, z]

def cart2sph(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    lat = np.arcsin(z/r)*180/np.pi
    lon = np.arctan2(y, x)*180/np.pi
    return [lat, lon, r]

def find_antipode(lat1, lon1, lat2, lon2):
    r = 6371.*1000.
    p1 = sph2cart(lat1, lon1, r)
    p2 = sph2cart(lat2, lon2, r)
    nrm = np.cross(p1, p2)

    new1 = cart2sph(nrm[0], nrm[1], nrm[2])
    new2 = cart2sph(-nrm[0], -nrm[1], -nrm[2])
    # return lat1, lon2 and lat2, lon2

    return new1[0], new1[1], new2[0], new2[1]







inv = client.get_stations(network='IU,IC,CU,II', station='*', starttime=eve, 
                endtime = eve + hours*60*60, level="response", location='00', channel=chan_code)


st = client.get_waveforms(network='IU,IC,CU,II', station='*', location='00', channel=chan_code, 
            starttime=eve, endtime=eve + hours*60*60)

st.merge(fill_value=0)
st.detrend('linear')
st.detrend('constant')
st.sort(['station'])


# Window
for tr in st:
    tr.data *= signal.get_window(('kaiser', 2.*np.pi), tr.stats.npts)

NFFT=2**(math.ceil(math.log(st[0].stats.npts, 2)))


lats, lons, amps, freqs, snrs = [], [], [], [], []
for tr in st:
    f,p = periodogram(tr.data, fs=tr.stats.sampling_rate, nfft= NFFT, scaling='spectrum')
    p, f = p[1:], f[1:]
    trid = tr.id
    inv_resp = inv.get_response(trid, tr.stats.starttime)
    resp, _ = inv_resp.get_evalresp_response(tr.stats.delta, NFFT, 'ACC')
    resp = resp[1:]
    # Convert units to nm/s/s

    p = np.sqrt(p/(np.abs(resp)**2))*10**9
    # Now have p in nm/s/s switch f to mHz
    f *= 1000.
    noise = p[(f>=3.0) & (f <=3.4)]
    p = p[(f >= mf) & (f <= Mf)]
    f = f[(f >= mf) & (f <= Mf)]
    
    snr = np.max(p)/np.max(noise)

    snrs.append(snr)
    snrs.append(snr)
    locs = inv.get_coordinates(tr.id, eve)
    lat1, lon1, lat2, lon2 = find_antipode(locs['latitude'], locs['longitude'], -20.5, 175.4)
    #blah = find_antipode(0, locs['longitude'], 0, 175.4)
    #print(blah)
    #import sys
    #sys.exit()

    lats.append(lat1)
    lats.append(lat2)
    lons.append(lon1)
    lons.append(lon2)
    amps.append(np.max(p))
    amps.append(np.max(p))
    freqs.append(f[np.argmax(p)])
    freqs.append(f[np.argmax(p)])



    #lats.append(locs['latitude'])
    #lons.append(locs['longitude'])
    #amps.append(np.max(p))
    #freqs.append(f[np.argmax(p)])
    if tr.stats.station in ['PFO', 'ADK']:
        print(tr.stats.station + ' ' + str(f[np.argmax(p)]))
print(lats)
print(lons)
gfres = []
gamps = []
for freq,snr  in zip(freqs,snrs):
    if snr >= 3:
        gfres.append(freq)
        gamps.append(amps)

fig = plt.figure(1, figsize=(14,14))
for idx  in range(2):
    ax = fig.add_subplot(1,2,1+idx)

    ax.scatter(175.4, -20.5, c='r',marker='*',s= 200,zorder=3, vmin=0., vmax=10)
    ax.scatter(180-175.4, 20.5, c='k',marker='*',s= 200, zorder=3, vmin=0., vmax=10)
    #ax.scatter(180-175.4, 20.5, c='k',marker='*',s= 200, transform=ccrs.Geodetic(), zorder=3, vmin=0., vmax=10)
    if idx ==0:
        var = amps
        label = 'Amplitude ($nm/s^2$)'
        lett = '(a)'
        minmax =[0,1]
    else:
        var = freqs
        label = 'Frequency (mHz)'
        lett = '(b)'
        minmax = [3.65,3.75]
    ax.set_title(lett, loc='left')
    for lat, lon, cvar, snr in zip(lats, lons, var, snrs):
        if snr >= 3:
            im = ax.scatter(lon, lat, c=cvar, s= 200, zorder=3, alpha=0.5,vmin=minmax[0],vmax=minmax[1], marker='o')
        else:
            im = ax.scatter(lon, lat, c=cvar, s= 50, zorder=3, alpha=0.5,vmin=minmax[0],vmax=minmax[1], marker='s')
    cbar = plt.colorbar(im, orientation='horizontal')
    cbar.set_label(label) 
    plt.xlabel('Longitude ($^{\circ}$)')
    plt.ylabel('Latitude ($^{\circ}$)')

plt.tight_layout()
plt.savefig('FigureS1.PNG', format='PNG', dpi=400)
plt.savefig('FigureS1.PDF', format='PDF', dpi=400)
plt.close('all') 



# print('Here is the frequency:' + str(np.mean(freqs)) + ' ' + str(np.std(freqs)))
# print('Here is the frequency with high SNR:' + str(np.mean(gfres)) + ' ' + str(np.std(gfres)))
# print('Here is the amps:' + str(np.mean(amps)) + ' ' + str(np.std(amps)))
# fig = plt.figure(1, figsize=(12,6))
# for idx  in range(2):
#     ax = fig.add_subplot(1,2,1+idx, projection = ccrs.Robinson())
#     ax.coastlines()
#     ax.add_feature(cart.feature.LAND, zorder=2, edgecolor='k')
#     ax.set_global()
#     ax.coastlines()
#     ax.scatter(175.4, -20.5, c='r',marker='*',s= 200, transform=ccrs.Geodetic(), zorder=3, vmin=0., vmax=10)
#     ax.scatter(180-175.4, 20.5, c='k',marker='*',s= 200, transform=ccrs.Geodetic(), zorder=3, vmin=0., vmax=10)
#     #ax.scatter(180-175.4, 20.5, c='k',marker='*',s= 200, transform=ccrs.Geodetic(), zorder=3, vmin=0., vmax=10)
#     if idx ==0:
#         var = amps
#         label = 'Amplitude ($nm/s^2$)'
#         lett = '(a)'
#         minmax =[0,1]
#     else:
#         var = freqs
#         label = 'Frequency (mHz)'
#         lett = '(b)'
#         minmax = [3.65,3.75]
#     ax.set_title(lett, loc='left')
#     for lat, lon, cvar, snr in zip(lats, lons, var, snrs):
#         if snr >= 3:
#             im = ax.scatter(lon, lat, c=cvar, s= 200, transform=ccrs.Geodetic(), zorder=3, alpha=0.5,vmin=minmax[0],vmax=minmax[1], marker='o')
#         else:
#             im = ax.scatter(lon, lat, c=cvar, s= 50, transform=ccrs.Geodetic(), zorder=3, alpha=0.5,vmin=minmax[0],vmax=minmax[1], marker='s')
#     cbar = plt.colorbar(im, orientation='horizontal')
#     cbar.set_label(label) 
# #plt.savefig('Tonga_' + chan + '.png', format='PNG', dpi=400)

# plt.tight_layout()
# plt.savefig('Figure5.PNG', format='PNG', dpi=400)
# plt.savefig('Figure5.PDF', format='PDF', dpi=400)
# plt.close('all') 

