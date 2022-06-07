from netCDF4 import Dataset
from datetime import *
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from glob import glob
from scipy.optimize import curve_fit
from mpl_toolkits.basemap import Basemap


'''

Calcula o e-folding scale dos dados altimetricos

'''


def func(t,a,b,c):
    return a*np.exp(-b*t) + c
lag = np.arange(200)

# path='/home/igor/Documents/AntSimi/altdata/'
# files=glob(path + '*.nc')
# files.sort()


# var=Dataset(files[ox])


plt.rcParams['font.size']=13
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.weight']='bold'

var=Dataset('/home/igor/Documents/AntSimi/altdata/clim/sla_climatology.nc')
lonc,latc=var['longitude'][:]-360,var['latitude'][:];
slac=var['sla'][:]
# adtc=var['adt'][:]

lat=var['latitude'][:]
lon=var['longitude'][:]-360

dt0=datetime(1950,1,1)
dt=[dt0 + timedelta(days=np.int(dia)) for dia in var['time']]

'E folding Scale Field'
EF=np.zeros([slac.shape[1],slac.shape[2]])



for iy in range(slac.shape[1]):
    for ix in range(slac.shape[2]):
        dsla=slac[:,iy,ix].data
        dsla=pd.Series(dsla)

        try:
            'Lagged Correlation'
            corr = np.hstack([dsla.corr(dsla.shift(i)) for i in lag])
            popt, pcov = curve_fit(func, lag, corr)
            y=func(lag,popt[0],popt[1],popt[2])


            EF[iy,ix]=1/popt[1]
        except:
            EF[iy,ix]=np.nan;





bvar=Dataset('/home/igor/Documents/Extras/Sumida/MERCATOR/bat/etopo1_bedrock.nc')
lonb=bvar['lon'][:];latb=bvar['lat'][:];bat=bvar['Band1'][:]
lonb,latb=np.meshgrid(lonb,latb)


'Mapa de E-Folding Scale'
lnmn=-52;lnmx=-31;ltmn=-29;ltmx=-10;

lonx,latx=np.meshgrid(lonc,latc)

p=plt.figure(figsize=(11.17,  8.77))
lbs=[True, False,False,True]
# C = Basemap(llcrnrlon=-49.5,llcrnrlat=-28,urcrnrlat=-23,urcrnrlon=-43.5,resolution='h',projection='merc')
C = Basemap(llcrnrlon=lnmn,llcrnrlat=ltmn,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);

C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,2),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn+2,ltmx,2),labels=lbs,zorder=101,linewidth=0.1)

cc=C.contourf(lonx,latx,EF,cmap='rainbow',levels=np.linspace(0,np.nanmax(EF),35),extend='max')

cbar=plt.colorbar(cc,format='%i',ticks=np.arange(0,280,40),aspect=50,pad=0.01)
cbar.set_label(r'E-folding rate [day]')
C.contourf(lonb,latb,bat,levels=[-150,np.inf],colors='dimgray',alpha=1,latlon=True)
C.contour(lonb,latb,-bat,levels=[-np.inf,150],colors='k',alpha=1,latlon=True,linewidhts=0.5)
ctt=C.contour(lonb,latb,-bat,levels=[1e3,2e3,3e3],colors='dimgray',alpha=0.8,latlon=True,linewidths=0.5)

plt.clabel(ctt,fmt='%i',fontsize=10)
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/clim/map_efold.png')
plt.close()





'Latitude'

p=plt.figure(figsize=(6.4 , 9.16))
plt.plot(np.nanmean(EF,axis=1),lat,'m')
plt.plot([20,200],[-20.5,-20.5],linestyle='--',color='orange')
plt.text(25,-20.4,'VTR',color='orange')
plt.xlim(20,200)
plt.ylabel('Latitude')
plt.xlabel('E-folding scale')
plt.grid()
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/lat_efold.png')




'Fit única'

Corr=np.zeros([slac.shape[1],slac.shape[2],lag.shape[0]])

for iy in range(slac.shape[1]):
    for ix in range(slac.shape[2]):
        dsla=slac[:,iy,ix].data
        dsla=pd.Series(dsla)
        corr = np.hstack([dsla.corr(dsla.shift(i)) for i in lag])
        Corr[iy,ix,:]=corr

Corr1=Corr.reshape(slac.shape[1]*slac.shape[2],lag.shape[0])


corrf=np.nanmean(Corr1,axis=0)
popt, pcov = curve_fit(func, lag, corrf)
yf=func(lag,popt[0],popt[1],popt[2])
plt.figure()
plt.plot(lag,corrf,'m',label='Correlation')
plt.plot(lag,yf,'c',label='Exponential Fit')
plt.title('e folding scale = {0:.2f}'.format(1/popt[1]))
plt.grid(alpha=0.5)


# 'INCOMPLETO'
# 'Latitude'

# Corr2=Corr.reshape(slac.shape[1],slac.shape[2],lag.shape[0])

# corrlat=np.nanmean(Corr2,axis=0)
# EFOLD=[]

# for ie in range(corrlat.shape[0]):
#     cl=corrlat[ie,:]
#     popt, pcov = curve_fit(func, lag, cl)
#     yf=func(lag,popt[0],popt[1],popt[2])
#     EFOLD.append(1/popt[1])

# plt.figure()
# plt.plot(EFOLD,lat,'m')
# plt.grid()

