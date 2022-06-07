from netCDF4 import Dataset
from datetime import *
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from glob import glob
from mpl_toolkits.basemap import Basemap
import seawater as sw
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy.ndimage.filters import gaussian_filter as gfilter
import sys
sys.path.append('/home/igor/Dropbox/Mestrado/Scripts/master/iury_wg/');
from scipy.interpolate import griddata
import gc as gc
from matplotlib.path import Path
import pandas as pd

'Bathymetry'
bvar=Dataset('/home/igor/Documents/Extras/Sumida/MERCATOR/bat/etopo1_bedrock.nc')
lonb=bvar['lon'][:];latb=bvar['lat'][:];bat=bvar['Band1'][:]
lonb,latb=np.meshgrid(lonb,latb)

'SLA and ADT'
var=Dataset('/home/igor/Documents/AntSimi/altdata/clim/adt_climatology.nc')
avar=Dataset('/home/igor/Documents/AntSimi/altdata/clim/sla_climatology.nc')

'Time'
time=[datetime(1950,1,1) + timedelta(days=np.int(dia)) for dia in var['time'][:].data]
time=pd.DatetimeIndex(time);
mo=time.month

'lat e lon'
lat=var['latitude'][:]
lon=var['longitude'][:]-360
lonx,latx=np.meshgrid(lon,lat)


'Delimited Area'
# dd=pd.read_csv('/home/igor/Documents/AntSimi/data/buffer_200m_6rad.csv')
dd=pd.read_csv('/home/igor/Documents/AntSimi/data/buffer_200m_200km.csv')
ddlon,ddlat=np.array(dd['lon']),np.array(dd['lat'])
ddseg=Path(np.array([ddlon,ddlat]).T)

maskx=np.empty(lonx.shape)*np.nan;

for yi in range(lonx.shape[1]):
    for xi in range(lonx.shape[0]):
        if ddseg.contains_point((lonx[yi,xi],latx[yi,xi])):
            maskx[yi,xi]=True


'Adt'
adtmask=var['adt'][:].mask
adt=var['adt'][:].data
adt[adtmask]=np.nan;

'mdt'
mdt=np.nanmean(adt,axis=0)
del adtmask

'Total velocity'
dx,dy=np.diff(lon).mean()*111e3,np.diff(lat).mean()*111e3
ut=-(9.8/sw.f(lat.data))*np.gradient(adt,axis=1)/dy
vt=(9.8/sw.f(lat.data))*np.gradient(adt,axis=2)/dx
del adt


'Sla'
slamask=avar['sla'][:].mask
sla=avar['sla'][:].data
sla[slamask]=np.nan;
del slamask

'Pertubation velocities'
up=-(9.8/sw.f(lat.data))*np.gradient(sla,axis=1)/dy
vp=(9.8/sw.f(lat.data))*np.gradient(sla,axis=2)/dx
del sla

'mdt velocity'
um=-(9.8/sw.f(lat.data))*np.gradient(mdt,axis=0)/dy
vm=(9.8/sw.f(lat.data))*np.gradient(mdt,axis=1)/dx

'Kinectic energies'
'MKE'
mke = 0.5*(um**2 + vm**2)

'EKE'
meke= 0.5*(up**2 + vp**2)
meke=np.nanmean(meke,axis=0)






'TOTAL MKE %%%%%%%%%%%%%%%%'

plt.rcParams['font.size']=20
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.weight']='bold'
ltmn=-29;ltmx=-10;lnmn=-52;lnmx=-35.

p=plt.figure(figsize=(6.54, 6.24))
lbs=[True, False,False,True]
# C = Basemap(llcrnrlon=-49.5,llcrnrlat=-28,urcrnrlat=-23,urcrnrlon=-43.5,resolution='h',projection='merc')
C = Basemap(llcrnrlon=lnmn+2,llcrnrlat=ltmn,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);

C.fillcontinents(color='#fecf78',zorder=100)
C.drawstates(linewidth=0.5,zorder=101)
C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,zorder=102,linewidth=0.1)
C.drawparallels(np.arange(ltmn+4,ltmx,4),labels=lbs,zorder=103,linewidth=0.1)

kratio1=mke
kratio1[maskx!=1]=np.nan;


cc=C.contourf(lonx,latx,kratio1,cmap='RdYlBu_r',levels=np.linspace(0,0.07,20),extend='max')
# C.contour(lonx,latx,((np.nanmean(meke,axis=0))/(np.nanmean(mke,axis=0)))*100,colors='navy',zorder=12,linewidths=1,levels=np.linspace(10,100,10))

# Q=C.quiver(lon,lat,und,vnd,scale=10)
# plt.quiverkey(Q,0.11,0.9,0.2,r'0.2 m s $^{-1}$')

cbar=plt.colorbar(cc,ticks=np.arange(0,0.08,0.01),format='%.2f',aspect=50,pad=0.02)
cbar.set_label(r'MKE [m$^{2}$s$^{-2}$]')
# # cbar.set_label(r'$\psi$ x 10$^{-5}$ [m$^2$ s$^{-1}$]')


C.contourf(lonb,latb,bat,levels=[-150,np.inf],colors='dimgray',alpha=1,latlon=True,zorder=14)
C.contour(lonb,latb,-bat,levels=[150],colors='k',alpha=1,linewidth=1,latlon=True,zorder=15)
# C.contour(lonb,latb,-bat,levels=[-np.inf,150],colors='k',alpha=1,latlon=True,linewidhts=0.5,zoder=16)
# ctt=C.contour(lonb,latb,-bat,levels=[1e3,2e3,3e3],colors='dimgray',alpha=0.8,latlon=True,linewidths=0.5)

# plt.clabel(ctt,fmt='%i',fontsize=10)
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/clim/MKE_total.png')
# plt.savefig('/home/igor/Documents/AntSimi/figures/clim/adt_clim.png')
plt.close()





'TOTAL mEKE %%%%%%%%%%%%%%%%'

plt.rcParams['font.size']=20
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.weight']='bold'

p=plt.figure(figsize=(6.54, 6.24))
lbs=[True, False,False,True]
# C = Basemap(llcrnrlon=-49.5,llcrnrlat=-28,urcrnrlat=-23,urcrnrlon=-43.5,resolution='h',projection='merc')
C = Basemap(llcrnrlon=lnmn+2,llcrnrlat=ltmn,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);

C.fillcontinents(color='#fecf78',zorder=100)
C.drawstates(linewidth=0.5,zorder=101)
C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,zorder=102,linewidth=0.1)
C.drawparallels(np.arange(ltmn+4,ltmx,4),labels=lbs,zorder=103,linewidth=0.1)

kratio2=meke
kratio2[maskx!=1]=np.nan;


cc=C.contourf(lonx,latx,kratio2,cmap='RdYlBu_r',levels=np.linspace(0,0.04,20),extend='max')
# C.contour(lonx,latx,((np.nanmean(meke,axis=0))/(np.nanmean(mke,axis=0)))*100,colors='navy',zorder=12,linewidths=1,levels=np.linspace(10,100,10))

# Q=C.quiver(lon,lat,und,vnd,scale=10)
# plt.quiverkey(Q,0.11,0.9,0.2,r'0.2 m s $^{-1}$')

cbar=plt.colorbar(cc,ticks=np.arange(0,0.06,0.01),format='%.2f',aspect=50,pad=0.02)
cbar.set_label(r'MEKE [m$^{2}$s$^{-2}$]')
# # cbar.set_label(r'$\psi$ x 10$^{-5}$ [m$^2$ s$^{-1}$]')


C.contourf(lonb,latb,bat,levels=[-150,np.inf],colors='dimgray',alpha=1,latlon=True,zorder=14)
C.contour(lonb,latb,-bat,levels=[150],colors='k',alpha=1,linewidth=1,latlon=True,zorder=15)
# C.contour(lonb,latb,-bat,levels=[-np.inf,150],colors='k',alpha=1,latlon=True,linewidhts=0.5,zoder=16)
# ctt=C.contour(lonb,latb,-bat,levels=[1e3,2e3,3e3],colors='dimgray',alpha=0.8,latlon=True,linewidths=0.5)

# plt.clabel(ctt,fmt='%i',fontsize=10)
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/clim/MEKE_total.png')
# plt.savefig('/home/igor/Documents/AntSimi/figures/clim/adt_clim.png')
plt.close()



print('North MEKE/MKE  = {0}'.format((np.nanmean(kratio2[45:,:])/np.nanmean(kratio1[45:,:]))*100))
print('South MEKE/MKE  = {0}'.format((np.nanmean(kratio2[:46,:])/np.nanmean(kratio1[:46,:]))*100))
print('Total  MEKE/MKE  = {0}'.format((np.nanmean(kratio2)/np.nanmean(kratio1))*100))


adtmask=var['adt'][:].mask
adt=var['adt'][:].data
adt[adtmask]=np.nan;

for ix in range(4):

    imo=[[1,2,3],
        [4,5,6],
        [7,8,9],
        [10,11,12]]
    seastr=['Summer','Autumn','Winter','Spring']


    cond=((mo>=imo[ix][0])&(mo<=imo[ix][-1]))

    sadt=np.nanmean(adt[cond,:,:],axis=0)
    ut=-(9.8/sw.f(lat.data))*np.gradient(sadt,axis=0)/dy
    vt=(9.8/sw.f(lat.data))*np.gradient(sadt,axis=1)/dx

    smke=0.5*(ut**2 + vt**2)

    kratio=smke


    'SEASONAL MKE'


    p=plt.figure(figsize=(6.54, 6.24))
    lbs=[True, False,False,True]
    # C = Basemap(llcrnrlon=-49.5,llcrnrlat=-28,urcrnrlat=-23,urcrnrlon=-43.5,resolution='h',projection='merc')
    C = Basemap(llcrnrlon=lnmn+2,llcrnrlat=ltmn,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

    C.drawcoastlines(linewidth=1);

    C.fillcontinents(color='#fecf78',zorder=100)
    C.drawstates(linewidth=0.5,zorder=101)
    C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,zorder=102,linewidth=0.1)
    C.drawparallels(np.arange(ltmn+4,ltmx,4),labels=lbs,zorder=103,linewidth=0.1)

    kratio1=smke
    kratio1[maskx!=1]=np.nan;


    cc=C.contourf(lonx,latx,kratio1,cmap='RdYlBu_r',levels=np.linspace(0,0.05,20),extend='max')
    # C.contour(lonx,latx,((np.nanmean(meke,axis=0))/(np.nanmean(mke,axis=0)))*100,colors='navy',zorder=12,linewidths=1,levels=np.linspace(10,100,10))

    # Q=C.quiver(lon,lat,und,vnd,scale=10)
    # plt.quiverkey(Q,0.11,0.9,0.2,r'0.2 m s $^{-1}$')

    cbar=plt.colorbar(cc,ticks=np.arange(0,0.08,0.01),format='%.2f',aspect=50,pad=0.02)
    cbar.set_label(r'MKE [m$^{2}$s$^{-2}$]')
    # # cbar.set_label(r'$\psi$ x 10$^{-5}$ [m$^2$ s$^{-1}$]')
    plt.title('{0}'.format(seastr[ix]))

    C.contourf(lonb,latb,bat,levels=[-150,np.inf],colors='dimgray',alpha=1,latlon=True,zorder=14)
    C.contour(lonb,latb,-bat,levels=[150],colors='k',alpha=1,linewidth=1,latlon=True,zorder=15)
    # C.contour(lonb,latb,-bat,levels=[-np.inf,150],colors='k',alpha=1,latlon=True,linewidhts=0.5,zoder=16)
    # ctt=C.contour(lonb,latb,-bat,levels=[1e3,2e3,3e3],colors='dimgray',alpha=0.8,latlon=True,linewidths=0.5)

    # plt.clabel(ctt,fmt='%i',fontsize=10)
    plt.tight_layout()
    plt.savefig('/home/igor/Documents/AntSimi/figures/clim/MKE_{0}.png'.format(seastr[ix]))
    # plt.savefig('/home/igor/Documents/AntSimi/figures/clim/adt_clim.png')
    plt.close()



    'SEASONAL MEKE'
    cond=((mo>=imo[ix][0])&(mo<=imo[ix][-1]))
    meke=0.5*(up**2 + vp**2)
    kratio2=np.nanmean(meke[cond,:,:],axis=0)
    kratio2[maskx!=1]=np.nan;




    p=plt.figure(figsize=(6.54, 6.24))
    lbs=[True, False,False,True]
    # C = Basemap(llcrnrlon=-49.5,llcrnrlat=-28,urcrnrlat=-23,urcrnrlon=-43.5,resolution='h',projection='merc')
    C = Basemap(llcrnrlon=lnmn+2,llcrnrlat=ltmn,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

    C.drawcoastlines(linewidth=1);

    C.fillcontinents(color='#fecf78',zorder=100)
    C.drawstates(linewidth=0.5,zorder=101)
    C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,zorder=102,linewidth=0.1)
    C.drawparallels(np.arange(ltmn+4,ltmx,4),labels=lbs,zorder=103,linewidth=0.1)




    cc=C.contourf(lonx,latx,kratio2,cmap='RdYlBu_r',levels=np.linspace(0,0.05,20),extend='max')
    # C.contour(lonx,latx,((np.nanmean(meke,axis=0))/(np.nanmean(mke,axis=0)))*100,colors='navy',zorder=12,linewidths=1,levels=np.linspace(10,100,10))

    # Q=C.quiver(lon,lat,und,vnd,scale=10)
    # plt.quiverkey(Q,0.11,0.9,0.2,r'0.2 m s $^{-1}$')

    cbar=plt.colorbar(cc,ticks=np.arange(0,0.08,0.01),format='%.2f',aspect=50,pad=0.02)
    cbar.set_label(r'MEKE [m$^{2}$s$^{-2}$]')
    # # cbar.set_label(r'$\psi$ x 10$^{-5}$ [m$^2$ s$^{-1}$]')


    C.contourf(lonb,latb,bat,levels=[-150,np.inf],colors='dimgray',alpha=1,latlon=True,zorder=14)
    C.contour(lonb,latb,-bat,levels=[150],colors='k',alpha=1,linewidth=1,latlon=True,zorder=15)
    # C.contour(lonb,latb,-bat,levels=[-np.inf,150],colors='k',alpha=1,latlon=True,linewidhts=0.5,zoder=16)
    # ctt=C.contour(lonb,latb,-bat,levels=[1e3,2e3,3e3],colors='dimgray',alpha=0.8,latlon=True,linewidths=0.5)
    plt.title('{0}'.format(seastr[ix]))

    # plt.clabel(ctt,fmt='%i',fontsize=10)
    plt.tight_layout()
    plt.savefig('/home/igor/Documents/AntSimi/figures/clim/MEKE_{0}.png'.format(seastr[ix]))
    # plt.savefig('/home/igor/Documents/AntSimi/figures/clim/adt_clim.png')
    plt.close()

    print('North MEKE/MKE {0} = {1}'.format(seastr[ix],(np.nanmean(kratio2[45:,:])/np.nanmean(kratio1[45:,:]))*100))
    print('South MEKE/MKE {0} = {1}'.format(seastr[ix],(np.nanmean(kratio2[:46,:])/np.nanmean(kratio1[:46,:]))*100))
    print('Total  MEKE/MKE {0} = {1}'.format(seastr[ix],(np.nanmean(kratio2)/np.nanmean(kratio1))*100))

del adt



