from netCDF4 import Dataset
from datetime import *
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from glob import glob
from mpl_toolkits.basemap import Basemap
import seawater as sw
import matplotlib.gridspec as gridspec
import scipy.stats as stats
from scipy.linalg import *
import seaborn as sns
from scipy.ndimage.filters import gaussian_filter as gfilter
import sys
sys.path.append('/home/igor/Dropbox/Mestrado/Scripts/master/iury_wg/');
from functions.lietal import uv2psiphi
from scipy.interpolate import griddata


'Plotagem de climatologias'



bvar=Dataset('/home/igor/Documents/Extras/Sumida/MERCATOR/bat/etopo1_bedrock.nc')
lonb=bvar['lon'][:];latb=bvar['lat'][:];bat=bvar['Band1'][:]
lonb,latb=np.meshgrid(lonb,latb)



var=Dataset('/home/igor/Documents/AntSimi/altdata/clim/adt_climatology.nc')
avar=Dataset('/home/igor/Documents/AntSimi/altdata/clim/sla_climatology.nc')







time=[datetime(1950,1,1) + timedelta(days=np.int(dia)) for dia in var['time'][:].data]
time=pd.DatetimeIndex(time);
mo=time.month

# right cut for the area
amask=var['adt'][:,11:,:].mask
adt=var['adt'][:,11:,:].data
adt[amask]=np.nan;
lat=var['latitude'][11:]
lon=var['longitude'][:]-360
dx,dy=np.diff(lon).mean()*111e3,np.diff(lat).mean()*111e3
lon,lat=np.meshgrid(lon,lat)

slamask=avar['sla'][:,11:,:].mask
sla=avar['sla'][:,11:,:].data
sla[slamask]=np.nan;



mdt=adt-sla
mdt=np.nanmean(mdt,axis=0)

'FIGURA 0: MDT'


plt.rcParams['font.size']=20
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.weight']='bold'
ltmn=-29;ltmx=-10;lnmn=-52;lnmx=-32.

p=plt.figure(figsize=(11.17,  8.77))
lbs=[True, False,False,True]
# C = Basemap(llcrnrlon=-49.5,llcrnrlat=-28,urcrnrlat=-23,urcrnrlon=-43.5,resolution='h',projection='merc')
C = Basemap(llcrnrlon=lnmn+2,llcrnrlat=ltmn,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);

C.fillcontinents(color='#fecf78',zorder=100)
C.drawstates(linewidth=0.5,zorder=101)
C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,zorder=102,linewidth=0.1)
C.drawparallels(np.arange(ltmn+4,ltmx,4),labels=lbs,zorder=103,linewidth=0.1)

# cc=C.contourf(lon,lat,-(9.8/sw.f(lat))*(hm-np.nanmean(hm))*1e-5,cmap='RdYlBu',levels=np.linspace(-2.5,2.5,35),zorder=10,extend='both')
# C.contour(lon,lat,(-9.8/sw.f(lat))*(hm-np.nanmean(hm))*1e-5,colors='navy',zorder=12,linewidths=1,levels=np.linspace(-2,2,35))

cc=C.contourf(lon,lat,mdt*100,cmap='rainbow',levels=np.linspace(40,70,35),extend='both')
ctt=C.contour(lon,lat,mdt*100,colors='navy',linewidths=1,levels=np.linspace(40,70,45))
cl=plt.clabel(ctt,fmt='%i',fontsize=10)
cl.zorder=13
# Q=C.quiver(lon,lat,und,vnd,scale=10)
# plt.quiverkey(Q,0.11,0.9,0.2,r'0.2 m s $^{-1}$')

cbar=plt.colorbar(cc,ticks=np.arange(-40,80,5),format='%.2f',aspect=50,pad=0.02)
cbar.set_label(r'$mdt$ [cm]')


# C.contourf(lonb,latb,bat,levels=[-150,np.inf],colors='dimgray',alpha=1,latlon=True,zorder=14)
# C.contour(lonb,latb,-bat,levels=[150],colors='k',alpha=1,linewidth=1,latlon=True,zorder=15)
# C.contour(lonb,latb,-bat,levels=[-np.inf,150],colors='k',alpha=1,latlon=True,linewidhts=0.5,zoder=16)
# ctt=C.contour(lonb,latb,-bat,levels=[1e3,2e3,3e3],colors='dimgray',alpha=0.8,latlon=True,linewidths=0.5)

plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/clim/mdt_clim.png')
plt.close()








hm=np.nanmean(adt,axis=0)



ux=-(9.8/sw.f(lat.data))*np.gradient(hm,axis=0)/dy
vx=(9.8/sw.f(lat.data))*np.gradient(hm,axis=1)/dx
u=ux.copy();v=vx.copy()


lonbx,latbx,batx=lonb[::10,::10],latb[::10,::10],bat[::10,::10]

meth='linear'
ui=griddata((lon[~np.isnan(ux)],lat[~np.isnan(ux)]),ux[~np.isnan(ux)],(lonbx,latbx),method=meth);
vi=griddata((lon[~np.isnan(vx)],lat[~np.isnan(vx)]),vx[~np.isnan(vx)],(lonbx,latbx),method=meth);

ui[np.isnan(ui)]=0;
vi[np.isnan(vi)]=0;

ui[batx>-200]=0;
vi[batx>-200]=0;



psix,und,vnd,_,_,_=uv2psiphi(lonbx,latbx,ui,vi)
psix=psix-np.nanmean(psix)

psix[np.isnan(ui)]=np.nan;
und[np.isnan(ui)]=np.nan;
vnd[np.isnan(vi)]=np.nan;






'FIGURA 1: CLIM V GEOST. PSI'

plt.rcParams['font.size']=20
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.weight']='bold'
ltmn=-29;ltmx=-10;lnmn=-52;lnmx=-32.

p=plt.figure(figsize=(11.17,  8.77))
lbs=[True, False,False,True]
# C = Basemap(llcrnrlon=-49.5,llcrnrlat=-28,urcrnrlat=-23,urcrnrlon=-43.5,resolution='h',projection='merc')
C = Basemap(llcrnrlon=lnmn+2,llcrnrlat=ltmn,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);

C.fillcontinents(color='#fecf78',zorder=100)
C.drawstates(linewidth=0.5,zorder=101)
C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,zorder=102,linewidth=0.1)
C.drawparallels(np.arange(ltmn+4,ltmx,4),labels=lbs,zorder=103,linewidth=0.1)

# cc=C.contourf(lon,lat,psix*1e-4,cmap='RdYlBu',levels=np.linspace(-2.5,2.5,25),zorder=10,extend='both')
# C.contour(lon,lat,psix*1e-4,colors='navy',zorder=12,linewidths=1,levels=np.linspace(-2,2,19))

cc=C.contourf(lonbx,latbx,psix*1e-4,cmap='RdYlBu',levels=np.linspace(-2.5,2.5,25),zorder=10,extend='both')
C.contour(lonbx,latbx,psix*1e-4,colors='navy',zorder=12,linewidths=1,levels=np.linspace(-2,2,19))

# Q=C.quiver(lon,lat,und,vnd,scale=10)
# plt.quiverkey(Q,0.11,0.9,0.2,r'0.2 m s $^{-1}$')

cbar=plt.colorbar(cc,ticks=np.arange(-2,3),format='%.2f',aspect=50,pad=0.02)
cbar.set_label(r'$\psi$ x 10$^{-4}$ [m$^2$ s$^{-1}$]')
# cbar.set_label(r'$\psi$ x 10$^{-5}$ [m$^2$ s$^{-1}$]')


C.contourf(lonb,latb,bat,levels=[-150,np.inf],colors='dimgray',alpha=1,latlon=True,zorder=14)
C.contour(lonb,latb,-bat,levels=[150],colors='k',alpha=1,linewidth=1,latlon=True,zorder=15)
C.contour(lonb,latb,-bat,levels=[-np.inf,150],colors='k',alpha=1,latlon=True,linewidhts=0.5,zoder=16)
# ctt=C.contour(lonb,latb,-bat,levels=[1e3,2e3,3e3],colors='dimgray',alpha=0.8,latlon=True,linewidths=0.5)

# plt.clabel(ctt,fmt='%i',fontsize=10)
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/clim/adt_clim+CC.png')
# plt.savefig('/home/igor/Documents/AntSimi/figures/clim/adt_clim.png')
plt.close()



'FIGURA 1b: CLIM V GEOST. VEL'



p=plt.figure(figsize=(11.17,  8.77))
lbs=[True, False,False,True]
# C = Basemap(llcrnrlon=-49.5,llcrnrlat=-28,urcrnrlat=-23,urcrnrlon=-43.5,resolution='h',projection='merc')
C = Basemap(llcrnrlon=lnmn+2,llcrnrlat=ltmn,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);

C.fillcontinents(color='#fecf78',zorder=100)
C.drawstates(linewidth=0.5,zorder=101)
C.drawmeridians(np.arange(lnmn,lnmx,2),labels=lbs,zorder=102,linewidth=0.1)
C.drawparallels(np.arange(ltmn+2,ltmx,2),labels=lbs,zorder=103,linewidth=0.1)

cc=C.contourf(lon,lat,np.sqrt(ux**2 + vx**2),cmap='rainbow',levels=np.linspace(0,0.35,35),zorder=10,extend='both')
# C.contour(lon,lat,np.sqrt(ux**2 + vx**2),colors='navy',zorder=12,linewidths=1,levels=np.linspace(-0.2,0.2,35))
dec=1
Q=C.quiver(lon[::dec,::dec],lat[::dec,::dec],ux[::dec,::dec],vx[::dec,::dec],scale=10,width=0.003,zorder=11)
plt.quiverkey(Q,0.11,0.9,0.2,r'0.2 m s $^{-1}$')

cbar=plt.colorbar(cc,ticks=np.arange(0,0.4,0.05),format='%.2f',aspect=50,pad=0.01)
# cbar.set_label(r'$\ADT$ x 10$^{-4}$ [m$^2$ s$^{-1}$]')
# cbar.set_label(r'$\psi$ x 10$^{-5}$ [m$^2$ s$^{-1}$]')


C.contourf(lonb,latb,bat,levels=[-150,np.inf],colors='dimgray',alpha=1,latlon=True,zorder=14)
C.contour(lonb,latb,-bat,levels=[150],colors='k',alpha=1,linewidth=1,latlon=True,zorder=15)
C.contour(lonb,latb,-bat,levels=[-np.inf,150],colors='k',alpha=1,latlon=True,linewidhts=0.5,zoder=16)
# ctt=C.contour(lonb,latb,-bat,levels=[1e3,2e3,3e3],colors='dimgray',alpha=0.8,latlon=True,linewidths=0.5)

# plt.clabel(ctt,fmt='%i',fontsize=10)
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/clim/adt_v.png')
plt.close()
del adt




'FIGURA 2: CLIM SLA RMS'
amask=avar['sla'][:,11:,:].mask

sla=avar['sla'][:,11:,:].data
sla[amask]=np.nan;

csla=sla.copy();csla[sla>0]=np.nan;
asla=sla.copy();asla[sla<0]=np.nan;

ham=sla #- np.nanmean(np.nanmean(sla,axis=1),axis=1)[:,None,None]

rms=np.sqrt(np.nanmean(ham**2,axis=0))

renq=np.sqrt(np.nanmean(ham**2,axis=0))/hm

del sla


lnmn=-52;lnmx=-31;ltmn=-29;ltmx=-10;

p=plt.figure(figsize=(11.17,  8.77))
lbs=[True, False,False,True]
# C = Basemap(llcrnrlon=-49.5,llcrnrlat=-28,urcrnrlat=-23,urcrnrlon=-43.5,resolution='h',projection='merc')
C = Basemap(llcrnrlon=lnmn+2,llcrnrlat=ltmn,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);

C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,2),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn+2,ltmx,2),labels=lbs,zorder=101,linewidth=0.1)

# cc=C.contourf(lon,lat,rms,cmap='jet',levels=np.linspace(0.01,np.nanmax(rms),35),extend='both')
cc=C.contourf(lon,lat,renq,cmap='jet',levels=np.linspace(0.05,0.15,35),extend='both')


cbar=plt.colorbar(cc,format='%.2f',aspect=50,pad=0.01)
# cbar.set_label(r'$RMSE$ [m]')
cbar.set_label(r'$NRMSE$ ')



C.contourf(lonb,latb,bat,levels=[-150,np.inf],colors='dimgray',alpha=1,latlon=True)
C.contour(lonb,latb,-bat,levels=[-np.inf,150],colors='k',alpha=1,latlon=True,linewidhts=0.5)
ctt=C.contour(lonb,latb,-bat,levels=[1e3,2e3,3e3],colors='dimgray',alpha=0.8,latlon=True,linewidths=0.5)

plt.clabel(ctt,fmt='%i',fontsize=10)
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/clim/sla_NRMS.png')
# plt.savefig('/home/igor/Documents/AntSimi/figures/clim/sla_RMS.png')

plt.close()


'FIGURA 3: CLIM SLA RMS CICLONICO'
sla=avar['sla'][:,11:,:].data
sla[amask]=np.nan;
csla=sla.copy();csla[sla>0]=np.nan;
del sla

cham=csla #- np.nanmean(np.nanmean(csla,axis=1),axis=1)[:,None,None]


crms=np.sqrt(np.nanmean(cham**2,axis=0))
del csla,cham



p=plt.figure(figsize=(11.17,  8.77))
lbs=[True, False,False,True]
# C = Basemap(llcrnrlon=-49.5,llcrnrlat=-28,urcrnrlat=-23,urcrnrlon=-43.5,resolution='h',projection='merc')
C = Basemap(llcrnrlon=lnmn+2,llcrnrlat=ltmn,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);

C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,2),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn+2,ltmx,2),labels=lbs,zorder=101,linewidth=0.1)

cc=C.contourf(lon,lat,crms,cmap='jet',levels=np.linspace(0,np.nanmax(crms),35),extend='both')


cbar=plt.colorbar(cc,format='%.2f',aspect=50,pad=0.01)
cbar.set_label(r'$RMS$ [m]')
C.contourf(lonb,latb,bat,levels=[-150,np.inf],colors='dimgray',alpha=1,latlon=True)
C.contour(lonb,latb,-bat,levels=[-np.inf,150],colors='k',alpha=1,latlon=True,linewidhts=0.5)
ctt=C.contour(lonb,latb,-bat,levels=[1e3,2e3,3e3],colors='dimgray',alpha=0.8,latlon=True,linewidths=0.5)

plt.clabel(ctt,fmt='%i',fontsize=10)
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/clim/cyclonic_RMS.png')
plt.close()



'FIGURA 4: CLIM SLA RMS ANTICICLONICO'

sla=avar['sla'][:,11:,:].data
sla[amask]=np.nan;
asla=sla.copy();asla[sla<0]=np.nan;
del sla

aham=asla - np.nanmean(np.nanmean(asla,axis=1),axis=1)[:,None,None]
del asla
arms=np.sqrt(np.nanmean(aham**2,axis=0))




p=plt.figure(figsize=(11.17,  8.77))
lbs=[True, False,False,True]
# C = Basemap(llcrnrlon=-49.5,llcrnrlat=-28,urcrnrlat=-23,urcrnrlon=-43.5,resolution='h',projection='merc')
C = Basemap(llcrnrlon=lnmn+2,llcrnrlat=ltmn,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);

C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,2),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn+2,ltmx,2),labels=lbs,zorder=101,linewidth=0.1)

cc=C.contourf(lon,lat,arms,cmap='jet',levels=np.linspace(0,np.nanmax(arms),35),extend='both')


cbar=plt.colorbar(cc,format='%.2f',aspect=50,pad=0.01)
cbar.set_label(r'$RMS$ [m]')
C.contourf(lonb,latb,bat,levels=[-150,np.inf],colors='dimgray',alpha=1,latlon=True)
C.contour(lonb,latb,-bat,levels=[-np.inf,150],colors='k',alpha=1,latlon=True,linewidhts=0.5)
ctt=C.contour(lonb,latb,-bat,levels=[1e3,2e3,3e3],colors='dimgray',alpha=0.8,latlon=True,linewidths=0.5)

plt.clabel(ctt,fmt='%i',fontsize=10)
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/clim/anticyclonic_RMS.png')
plt.close()



'FIGURA 5: CLIM ADT SEASON'

amask=var['adt'][:,11:,:].mask
adt=var['adt'][:,11:,:].data
adt[amask]=np.nan;

vm=np.nanmean(adt[mo<=3,:,:],axis=0)
om=np.nanmean(adt[(mo>3)&(mo<=6),:,:],axis=0)
im=np.nanmean(adt[(mo>7)&(mo<=9),:,:],axis=0)
pm=np.nanmean(adt[mo>=10,:,:],axis=0)

season=np.array([vm,om,im,pm])
seasonstr=['Summer','Fall','Winter','Spring']

del adt

for ixx in range(4):

    hm=season[ixx]
    u=-(9.8/sw.f(lat.data))*np.gradient(hm,axis=0)/dy
    v=(9.8/sw.f(lat.data))*np.gradient(hm,axis=1)/dx

    u[np.isnan(u)]=0;
    v[np.isnan(v)]=0;


    psix,und,vnd,_,_,_=uv2psiphi(lon,lat,u,v)
    psix=psix-np.nanmean(psix)

    p=plt.figure(figsize=(11.17,  8.77))
    lbs=[True, False,False,True]
    # C = Basemap(llcrnrlon=-49.5,llcrnrlat=-28,urcrnrlat=-23,urcrnrlon=-43.5,resolution='h',projection='merc')
    C = Basemap(llcrnrlon=lnmn+2,llcrnrlat=ltmn,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

    C.drawcoastlines(linewidth=1);

    C.fillcontinents(color='#fecf78')
    C.drawstates(linewidth=0.5)
    C.drawmeridians(np.arange(lnmn,lnmx,2),labels=lbs,zorder=100,linewidth=0.1)
    C.drawparallels(np.arange(ltmn+2,ltmx,2),labels=lbs,zorder=101,linewidth=0.1)

    cc=C.contourf(lon,lat,psix*1e-4,cmap='RdYlBu',levels=np.linspace(-2.5,2.5,35),extend='both')
    Q=C.quiver(lon,lat,und,vnd,scale=10)
    plt.quiverkey(Q,0.11,0.9,0.2,r'0.2 m s $^{-1}$')

    cbar=plt.colorbar(cc,ticks=np.arange(-2,3),format='%.2f',aspect=50,pad=0.01)
    cbar.set_label(r'$\psi$ x 10$^{-4}$ [m$^2$ s$^{-1}$]')
    C.contourf(lonb,latb,bat,levels=[-150,np.inf],colors='dimgray',alpha=1,latlon=True)
    C.contour(lonb,latb,-bat,levels=[-np.inf,150],colors='k',alpha=1,latlon=True,linewidhts=0.5)
    ctt=C.contour(lonb,latb,-bat,levels=[1e3,2e3,3e3],colors='dimgray',alpha=0.8,latlon=True,linewidths=0.5)

    plt.title('{0} Climatology'.format(seasonstr[ixx]),fontweight='bold')

    plt.clabel(ctt,fmt='%i',fontsize=10)
    plt.tight_layout()
    plt.savefig('/home/igor/Documents/AntSimi/figures/clim/adt_clim_{0}.png'.format(seasonstr[ixx]))
    plt.close()
