#!/usr/bin/env python
import matplotlib.pyplot as plt
from scipy.signal import periodogram
import numpy as np
import math
from scipy import signal
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
import cartopy.crs as ccrs
import cartopy as cart
import matplotlib as mpl
mpl.rc('font', family='serif')
mpl.rc('font', serif='Times')
#mpl.rc('text', usetex=True)
mpl.rc('font', size=18)



# Now let's see if we can invert for Q

# network
client = Client("IRIS")
eve = UTCDateTime('2022-01-15T04:14:00') 
#eve = UTCDateTime('1991-06-15T06:30:00')
mf, Mf = 3.65, 3.78
hours = 30
window = 15
timestep = 1



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
                          endtime=eve + hours*60*60, level="response",
                          location='00', channel='LHZ')

Qsg, resisg, latsg, lonsg = [], [], [], []
Qsb, resisb, latsb, lonsb = [], [], [], []
for net in inv:
    for sta in net:
        try:
            st = client.get_waveforms(network='IU,IC,CU,II', station=sta.code, location='00',
                                  channel='LHZ', starttime=eve,
                                  endtime=eve + hours*60*60)
        except:
            continue
        st.merge(fill_value=0)
        st.detrend('linear')
        vals = []

        for st_temp in st.slide(window*60*60, timestep*60*60):
            tr = st_temp[0].copy()
            tr.data *= signal.get_window(('kaiser', 2. * np.pi), tr.stats.npts)
            NFFT = 2 ** (math.ceil(math.log(tr.stats.npts, 2)))
            f, p = periodogram(tr.data, fs=tr.stats.sampling_rate,
                               nfft=NFFT, scaling='spectrum')
            p, f = p[1:], f[1:]
            # Here we are just selecting values of the spectrum around the 3.72mz peak
            f *= 1000.
            p = p[(f >= mf) & (f <= Mf)]
            f = f[(f >= mf) & (f <= Mf)]
            val = (1./(2*np.pi*Mf/1000.- 2*np.pi*mf/1000.))*np.trapz(p, 2*np.pi*f/1000.)
            # Just take the ln of the value
            vals.append(np.log(val))


        times = np.arange(len(vals)) + window/2.

        # Rob calculates cycles per time step, so we don't need this part (L86)
        times *= 60*60
        # fit a curve
        p, resi,_,_,_  = np.polyfit(times, vals, 1, full=True)
        Q = (2*np.pi)/(np.abs(p[0])*(1./((Mf/1000.+mf/1000.)/2.)))
        goodp = 1./(191.*(1./((Mf/1000.+mf/1000.)/2.))/(2*np.pi))

        print(sta.code + ' Q=' + str(Q) + ' Qerror' + str(resi))
        if resi <=1:
            if Q <=600:
                Qsg.append(Q)
                Qsg.append(Q)
            else:
                Qsg.append(600.)
                Qsg.append(600.)
            resisg.append(resi)
            resisg.append(resi)
            lat1, lon1, lat2, lon2 = find_antipode(sta.latitude, sta.longitude, -20.5, 175.4)
            latsg.append(lat1)
            latsg.append(lat2)
            lonsg.append(lon1)
            lonsg.append(lon2)
        else:
            Qsb.append(Q)
            resisb.append(resi)
            Qsb.append(Q)
            resisb.append(resi)
            lat1, lon1, lat2, lon2 = find_antipode(sta.latitude, sta.longitude, -20.5, 175.4)
            latsb.append(lat1)
            lonsb.append(lon1)
            latsb.append(lat2)
            lonsb.append(lon2)
        p1 = np.poly1d(p)
        vals = np.array(vals)
        print('Here we are')
        print(p)
        newvals = -goodp*times 
        print(newvals)
        print(10**(vals))
        fig = plt.figure(2, figsize=(12,12))
        plt.semilogy(times/(60*60), np.exp(vals),'.', label='Measurements')
        plt.semilogy(times/(60*60), p1(times), label='Q=' + str(round(Q,2)))
        plt.semilogy(times/(60*60),newvals , label='Q=191 PREM')
        plt.xlabel('Time (hr.)')
        plt.ylabel('Average Power (counts${}}^2$)')
        plt.legend()
        plt.show()
        plt.clf()


fig = plt.figure(1, figsize=(12,12))
ax = fig.add_subplot(1,1,1)

ax.scatter(175.4, -20.5, c='r',marker='*',s= 200, zorder=3, vmin=0., vmax=10)
ax.scatter(180-175.4, 20.5, c='k',marker='*',s= 200, zorder=3, vmin=0., vmax=10)
    # if idx ==0:
    #     var = amps
    #     label = 'Amplitude ($nm/s^2$)'
    #     lett = '(a)'
    #     minmax =[0,1]
    # else:
    #     var = freqs
    #     label = 'Frequency (mHz)'
    #     lett = '(b)'
    #     minmax = [3.65,3.75]
#ax.set_title(lett, loc='left')

im = ax.scatter(lonsb, latsb, c=Qsb, s = 50., zorder=3, alpha=0.5, vmin=100, vmax=600, marker='s')
im = ax.scatter(lonsg, latsg, c=Qsg, s = 200., zorder=3, alpha=0.7, vmin=100, vmax=600, marker='o')
cbar = plt.colorbar(im, orientation='horizontal')
cbar.set_label('Estimated Q') 
#plt.savefig('Tonga_' + chan + '.png', format='PNG', dpi=400)
plt.xlabel('Longitude ($^{\circ}$)')
plt.ylabel('Latitude ($^{\circ}$)')
plt.tight_layout()
plt.savefig('FigureS4.PNG', format='PNG', dpi=400)
plt.savefig('FigureS4.PDF', format='PDF', dpi=400)
plt.close('all') 

print(str(np.mean(Qsg)) + ' ' + str(np.std(Qsg)) )


