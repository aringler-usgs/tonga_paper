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
            inv_resp = inv.get_response(tr.id, tr.stats.starttime)
            resp, _ = inv_resp.get_evalresp_response(tr.stats.delta, NFFT, 'ACC')
            resp = resp[1:]
            p /= np.abs(resp)
            p *= 10**9
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




        newp = -((2*np.pi)/(191.))/(1./((Mf/1000.+mf/1000.)/2.))
        print(Q)
        print(sta.code + ' Q=' + str(Q) + ' Qerror' + str(resi))
        if resi <=1:
            if Q <=600:
                Qsg.append(Q)
            else:
                Qsg.append(600.)
            resisg.append(resi)
            latsg.append(sta.latitude)
            lonsg.append(sta.longitude)
        else:
            Qsb.append(Q)
            resisb.append(resi)
            latsb.append(sta.latitude)
            lonsb.append(sta.longitude)
        print('HEre is p')
        p2 = p
        print(p)
        print('Here is newp')
        print(newp)
        p = np.poly1d(p)
        #print(p)
        vals = np.array(vals)
        fig = plt.figure(2, figsize=(12,12))
        plt.semilogy(times/(60*60), np.exp(vals),'.', label='Measurements')
        plt.semilogy(times/(60*60), np.exp(p(times)), label='Q=' + str(Q))
        print('BOOOOOOOM')
        print(np.exp(np.mean(vals)))
        plt.semilogy(times/(60*60),np.exp(newp*times - (np.mean(vals) + (newp*times)[0])), label='Q=191 PREM')
        plt.xlabel('Time (hr.)')
        plt.ylabel('Average Power ($nm/s^2$)')
        plt.legend()
        plt.savefig('FigureS2.png', format='PNG', dpi=400)
        plt.savefig('FigureS2.pdf', format='PDF', dpi=400)
        import sys
        sys.exit()
        plt.show()
        plt.clf()


fig = plt.figure(1, figsize=(12,6))
ax = fig.add_subplot(1,1,1, projection = ccrs.Robinson())
ax.coastlines()
ax.add_feature(cart.feature.LAND, zorder=2, edgecolor='k')
ax.set_global()
ax.coastlines()
ax.scatter(175.4, -20.5, c='r',marker='*',s= 200, transform=ccrs.Geodetic(), zorder=3, vmin=0., vmax=10)
ax.scatter(180-175.4, 20.5, c='k',marker='*',s= 200, transform=ccrs.Geodetic(), zorder=3, vmin=0., vmax=10)
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

im = ax.scatter(lonsb, latsb, c=Qsb, s = 50., transform=ccrs.Geodetic(), zorder=3, alpha=0.5, vmin=100, vmax=600, marker='s')
im = ax.scatter(lonsg, latsg, c=Qsg, s = 200., transform=ccrs.Geodetic(), zorder=3, alpha=0.7, vmin=100, vmax=600, marker='o')
cbar = plt.colorbar(im, orientation='horizontal')
cbar.set_label('Estimated Q') 
#plt.savefig('Tonga_' + chan + '.png', format='PNG', dpi=400)

plt.tight_layout()
plt.savefig('Figure6.PNG', format='PNG', dpi=400)
plt.savefig('Figure6.PDF', format='PDF', dpi=400)
plt.close('all') 

print(str(np.mean(Qsg)) + ' ' + str(np.std(Qsg)) )


