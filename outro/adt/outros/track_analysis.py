from netCDF4 import Dataset
from py_eddy_tracker.dataset.grid import RegularGridDataset, UnRegularGridDataset
# import logging
# logging.basicConfig(level=logging.DEBUG)
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
import cmocean as cm

plt.rcParams['font.size']=13
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.weight']='bold'

var=Dataset('/home/igor/Documents/AntSimi/altdata/clim/sla_climatology.nc')

lonc,latc=var['longitude'][:]-360,var['latitude'][:];
# slac=var['sla'][:]

dt0=datetime(1950,1,1)
dt=[dt0 + timedelta(days=np.int(dia)) for dia in var['time']]

'Topography'
topo=Dataset('/home/igor/Documents/Extras/Sumida/HYCOM/bat/etopo1_bedrock.nc')
xb=topo['lon'][:];yb=topo['lat'][:];zb=topo['Band1'][:];
xxb,yyb=np.meshgrid(xb,yb)



anc=Dataset('/home/igor/Documents/AntSimi/tracking/Anticyclonic.nc')
cnc=Dataset('/home/igor/Documents/AntSimi/tracking/Cyclonic.nc')



atrack=anc['track'][:]
alon,alat=anc['longitude'][:]-360,anc['latitude'][:]

ctrack=cnc['track'][:]
clon,clat=cnc['longitude'][:]-360,cnc['latitude'][:]


ltmn=-32;ltmx=-9;lnmn=-52;lnmx=-30.

p=plt.figure(figsize=(10.68,  9.22))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=ltmn,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,2),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,2),labels=lbs,zorder=101,linewidth=0.1)

# [C.plot(alon[atrack==iy],alat[atrack==iy],'r',linewidth=0.5) for iy in range(atrack.max())]
[C.plot(alon[atrack==iy][0],alat[atrack==iy][0],'ro',markersize=1) for iy in range(atrack.max())]

# [C.plot(clon[ctrack==iy],clat[ctrack==iy],'b',linewidth=0.5) for iy in range(ctrack.max())]
[C.plot(clon[ctrack==iy][0],clat[ctrack==iy][0],'bo',markersize=1) for iy in range(ctrack.max())]

plt.savefig('/home/igor/Documents/AntSimi/figures/formation_points.png')
plt.close()



'Snapshot Tracked'
ref=anc['time'][:][0]
at=np.int_((anc['time'][:]-ref))
ct=np.int_((cnc['time'][:]-ref))

'''
for ip in range(slac.shape[0]):
    p=plt.figure(figsize=(6.4 , 8.76))
    plt.contourf(lonc,latc,slac[ip],cmap='seismic',levels=np.linspace(-0.5,0.5,100))

    [plt.plot(xll,yll,'ro',markersize=3) for xll,yll in zip(alon[at==ip],alat[at==ip]) ]
    [plt.text(xll,yll,np.int(tll),color='r') for xll,yll,tll in zip(alon[at==ip],alat[at==ip],atrack[at==ip]) ]

    [plt.plot(xll,yll,'bo',markersize=3) for xll,yll in zip(clon[ct==ip],clat[ct==ip]) ]
    [plt.text(xll,yll,np.int(tll),color='b') for xll,yll,tll in zip(clon[ct==ip],clat[ct==ip],ctrack[ct==ip]) ]
    
    plt.contour(xxb,yyb,zb,colors='k',alpha=1,levels=[-100,np.inf])
    plt.ylim(latc.min(),latc.max())
    plt.xlim(lonc.min(),lonc.max())

    plt.tight_layout()
    plt.savefig('/home/igor/Documents/AntSimi/maps/{0:04d}_track.png'.format(ip+1))
    plt.close()
'''

'Anticyclonic parameters'
Atime=anc['time'][:]
Aamp=pd.Series(anc['amplitude'][:]);Aamp=Aamp.groupby(atrack).mean();
Avel=pd.Series(anc['speed_average'][:]);Avel=Avel.groupby(atrack).mean();
Arad=pd.Series(anc['effective_radius'][:]);Arad=Arad.groupby(atrack).mean();
Aeke=0.5*((Avel*100)**2)

Adist=[];Adur=[];Apvel=[];Apvel2=[];
Adate=[];Afdate=[];


for ai in range(atrack.max()+1):

    dur=len(atrack[atrack==ai])

    distx=sw.dist(lon=alon[atrack==ai],lat=alat[atrack==ai])[0]
    k=np.insert(distx,0,0)
    dist=np.cumsum(k)

    pvel=(dist[-1]*1e3)/(dur*86400)
    pvel2=np.nanmean((distx*1e3)/(86400))

    time=Atime[atrack==ai];

    date=dt0+ timedelta(days=np.int(time.mean()));
    fdate=dt0 + timedelta(days=np.int(time[0]));


    Adist.append(dist[-1]);
    Adur.append(dur);
    Apvel.append(pvel);
    Apvel2.append(pvel2);
    Adate.append(date);
    Afdate.append(fdate);

Adist=np.array(Adist);
Adur=np.array(Adur);
Apvel=np.array(Apvel);
Apvel2=np.array(Apvel2);


'Cyclonic parameters'
Ctime=cnc['time'][:]
Camp=pd.Series(cnc['amplitude'][:]);Camp=Camp.groupby(ctrack).mean();
Cvel=pd.Series(cnc['speed_average'][:]);Cvel=Cvel.groupby(ctrack).mean();
Crad=pd.Series(cnc['effective_radius'][:]);Crad=Crad.groupby(ctrack).mean();
Ceke=0.5*((Cvel*100)**2)

Cdist=[];Cdur=[];Cpvel=[];Cpvel2=[];
Cdate=[];Cfdate=[];

for ci in range(ctrack.max()+1):

    dur=len(ctrack[ctrack==ci])

    distx=sw.dist(lon=clon[ctrack==ci],lat=clat[ctrack==ci])[0]

    k=np.insert(distx,0,0)
    dist=np.cumsum(k)

    pvel=(dist[-1]*1e3)/(dur*86400)
    pvel2=np.nanmean((distx*1e3)/86400)

    time=Ctime[ctrack==ci];

    date=dt0+ timedelta(days=np.int(time.mean()));
    fdate=dt0 + timedelta(days=np.int(time[0]));

    Cdist.append(dist[-1]);
    Cdur.append(dur);
    Cpvel.append(pvel);
    Cpvel2.append(pvel2);
    Cdate.append(date);
    Cfdate.append(fdate);


Cdist=np.array(Cdist)
Cdur=np.array(Cdur)
Cpvel=np.array(Cpvel)
Cpvel2=np.array(Cpvel2)


vamp=1
vrad=40


ccond=(Camp>=vamp/100)&(Crad>=vrad*1e3)
'Table'
ctt=pd.Series({'amp':np.nanmean(Camp[ccond]),'amp_min':np.nanmin(Camp[ccond]),'amp_max':np.nanmax(Camp[ccond]),'amp_std':np.nanstd(Camp[ccond]),
                'vel':np.nanmean(Cvel[ccond]),'vel_min':np.nanmin(Cvel[ccond]),'vel_max':np.nanmax(Cvel[ccond]),'vel_std':np.nanstd(Cvel[ccond]),
                'rad':np.nanmean(Crad[ccond]),'rad_min':np.nanmin(Crad[ccond]),'rad_max':np.nanmax(Crad[ccond]),'rad_std':np.nanstd(Crad[ccond]),
                'eke':np.nanmean(Ceke[ccond]),'eke_min':np.nanmin(Ceke[ccond]),'eke_max':np.nanmax(Ceke[ccond]),'eke_std':np.nanstd(Ceke[ccond]),
                'dist':np.nanmean(Cdist[ccond]),'dist_min':np.nanmin(Cdist[ccond]),'dist_max':np.nanmax(Cdist[ccond]),'dist_std':np.nanstd(Cdist[ccond]),
                'dur':np.nanmean(Cdur[ccond]),'dur_min':np.nanmin(Cdur[ccond]),'dur_max':np.nanmax(Cdur[ccond]),'dur_std':np.nanstd(Cdur[ccond]),
                'pvel':np.nanmean(Cpvel[ccond]),'pvel_min':np.nanmin(Cpvel[ccond]),'pvel_max':np.nanmax(Cpvel[ccond]),'pvel_std':np.nanstd(Cpvel[ccond]),
                'pvel2':np.nanmean(Cpvel2[ccond]),'pvel2_min':np.nanmin(Cpvel2[ccond]),'pvel2_max':np.nanmax(Cpvel2[ccond]),'pvel2_std':np.nanstd(Cpvel2[ccond])})

acond=(Aamp>=vamp/100)&(Arad>=vrad*1e3)

att=pd.Series({'amp':np.nanmean(Aamp[acond]),'amp_min':np.nanmin(Aamp[acond]),'amp_max':np.nanmax(Aamp[acond]),'amp_std':np.nanstd(Aamp[acond]),
                'vel':np.nanmean(Avel[acond]),'vel_min':np.nanmin(Avel[acond]),'vel_max':np.nanmax(Avel[acond]),'vel_std':np.nanstd(Avel[acond]),
                'rad':np.nanmean(Arad[acond]),'rad_min':np.nanmin(Arad[acond]),'rad_max':np.nanmax(Arad[acond]),'rad_std':np.nanstd(Arad[acond]),
                'eke':np.nanmean(Aeke[acond]),'eke_min':np.nanmin(Aeke[acond]),'eke_max':np.nanmax(Aeke[acond]),'eke_std':np.nanstd(Aeke[acond]),
                'dist':np.nanmean(Adist[acond]),'dist_min':np.nanmin(Adist[acond]),'dist_max':np.nanmax(Adist[acond]),'dist_std':np.nanstd(Adist[acond]),
                'dur':np.nanmean(Adur[acond]),'dur_min':np.nanmin(Adur[acond]),'dur_max':np.nanmax(Adur[acond]),'dur_std':np.nanstd(Adur[acond]),
                'pvel':np.nanmean(Apvel[acond]),'pvel_min':np.nanmin(Apvel[acond]),'pvel_max':np.nanmax(Apvel[acond]),'pvel_std':np.nanstd(Apvel[acond]),
                'pvel2':np.nanmean(Apvel2[acond]),'pvel2_min':np.nanmin(Apvel2[acond]),'pvel2_max':np.nanmax(Apvel2[acond]),'pvel2_std':np.nanstd(Apvel2[acond])})



savep='/home/igor/Documents/AntSimi/figures/'
ctt.to_csv(savep+ 'cyclone_prop.dat',sep='\t',mode='w',float_format='%.6f',index=True,header=True,na_rep='NaN')
att.to_csv(savep + 'anticyclone_prop.dat',sep='\t',mode='w',float_format='%.6f',index=True,header=True,na_rep='NaN')




'''
Figure 1a & 1b - Regression Plots
CC
'''

Wc=Crad/1e3

colorx='navy'
fsz=14;

pc=plt.figure(figsize=(20.01,  4.76))

gs = gridspec.GridSpec(1,19)
ax1 = plt.subplot(gs[0, 0:4])

Xc,Yc=-Camp*100,Cvel*100
# X=Xc[Xc<=-1];Y=Yc[Xc<=-1];
X=Xc[(Xc<=-vamp)&(Wc>=vrad)];Y=Yc[(Xc<=-vamp)&(Wc>=vrad)];


slope, intercept, r_value, p_value, std_err = stats.linregress(X,Y)
line = slope*(X)+intercept

plt.scatter(X,Y,c=colorx,edgecolors='white',linewidths=1,alpha=0.7)
plt.plot(X,line,color=colorx,alpha=0.5)

textstr = '\n'.join(('pearson r ={:.2f}'.format(r_value),
                    'p value ={:.2E}'.format(p_value)))

plt.title(textstr,fontweight='bold')

plt.grid()

plt.ylabel(r'Swirling Speed (cm s$^{-1}$)',fontsize=fsz, fontweight='bold')
plt.xlabel('Amplitude (cm)',fontsize=fsz, fontweight='bold')
plt.ylim(0,60)
plt.xlim(-20,0)



ax2 = plt.subplot(gs[0, 5:9])

Xc,Yc=-Camp*100,Ceke
# X=Xc[Xc<=-1];Y=Yc[Xc<=-1];
X=Xc[(Xc<=-vamp)&(Wc>=vrad)];Y=Yc[(Xc<=-vamp)&(Wc>=vrad)];

slope, intercept, r_value, p_value, std_err = stats.linregress(X,Y)
line = slope*(X)+intercept

plt.scatter(X,Y,c=colorx,edgecolors='white',linewidths=1,alpha=0.7)
plt.plot(X,line,color=colorx,alpha=0.5)

textstr = '\n'.join(('pearson r ={:.2f}'.format(r_value),
                    'p value ={:.2E}'.format(p_value)))

plt.title(textstr,fontweight='bold')
plt.grid()

plt.ylabel(r'EKE (cm$^2$ s$^{-2}$)',fontsize=fsz, fontweight='bold')
plt.xlabel('Amplitude (cm)',fontsize=fsz, fontweight='bold')
plt.ylim(0,Y.max())
plt.xlim(-20,0)


ax3 = plt.subplot(gs[0, 10:14])

Xc,Yc=-Camp*100, Crad/1e3
# X=Xc[Xc<=-1];Y=Yc[Xc<=-1];
# X=Xc[(Xc<=-1 )&( Yc>=40)];Y=Yc[(Xc<=-1 )&(Yc>=40)];
X=Xc[(Xc<=-vamp)&(Wc>=vrad)];Y=Yc[(Xc<=-vamp)&(Wc>=vrad)];




slope, intercept, r_value, p_value, std_err = stats.linregress(X,Y)
line = slope*(X)+intercept

plt.scatter(X,Y,c=colorx,edgecolors='white',linewidths=1,alpha=0.7)
plt.plot(X,line,color=colorx,alpha=0.5)

textstr = '\n'.join(('pearson r ={:.2f}'.format(r_value),
                    'p value ={:.2E}'.format(p_value)))


# plt.text(X.min(),Y.max()-0.5,textstr)
plt.title(textstr,fontweight='bold')
plt.grid()

plt.ylabel(r'Radius (km)',fontsize=fsz, fontweight='bold')
plt.xlabel('Amplitude (cm)',fontsize=fsz, fontweight='bold')
plt.ylim(30,150)
plt.xlim(-20,0)


ax4 = plt.subplot(gs[0,15:19])

xxc=-Camp*100
Xc,Yc=-Crad/1e3, Cpvel2*100
# X=Xc[Xc<=-1];Y=Yc[Xc<=-1];
X=Xc[(xxc<=-vamp)&(Wc>=vrad)];Y=Yc[(xxc<=-vamp)&(Wc>=vrad)];



slope, intercept, r_value, p_value, std_err = stats.linregress(X,Y)
line = slope*(X)+intercept

plt.scatter(X,Y,c=colorx,edgecolors='white',linewidths=1,alpha=0.7)
plt.plot(X,line,color=colorx,alpha=0.5)

textstr = '\n'.join(('pearson r ={:.2f}'.format(r_value),
                    'p value ={:.2E}'.format(p_value)))


# plt.text(X.min(),Y.max()-0.5,textstr)
plt.title(textstr,fontweight='bold')
plt.grid()

plt.ylabel(r'Propagating Speed (cm $s^{-1}$)',fontsize=fsz, fontweight='bold')
plt.xlabel('Radius (km)',fontsize=fsz, fontweight='bold')
plt.ylim(0,30)
plt.xlim(-150,-30)

plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/Scatter_Cyclone.png')
plt.close()



'''
AC
'''
Wc=Arad/1e3
colorx='red'
fsz=14;
pa=plt.figure(figsize=(20.01,  4.76))


gs = gridspec.GridSpec(1,19)
ax1 = plt.subplot(gs[0, 0:4])

Xc,Yc=Aamp*100,Avel*100

# X=Xc[Xc>=1];Y=Yc[Xc>=1];
# X=Xc[Xc>=np.nanmean(Xc)-np.nanstd(Xc)];Y=Yc[Xc>=np.nanmean(Xc)-np.nanstd(Xc)];
X=Xc[(Xc>=vamp)&(Wc>=vrad)];Y=Yc[(Xc>=vamp)&(Wc>=vrad)];

slope, intercept, r_value, p_value, std_err = stats.linregress(X,Y)
line = slope*(X)+intercept

plt.scatter(X,Y,c=colorx,edgecolors='white',linewidths=1,alpha=0.7)
plt.plot(X,line,color=colorx,alpha=0.5)

textstr = '\n'.join(('pearson r ={:.2f}'.format(r_value),
                    'p value ={:.2E}'.format(p_value)))


# plt.text(X.min(),Y.max()-0.5,textstr)
plt.title(textstr,fontweight='bold')
plt.grid()

plt.ylabel(r'Swirling Speed (cm s$^{-1}$)',fontsize=fsz, fontweight='bold')
plt.xlabel('Amplitude (cm)',fontsize=fsz, fontweight='bold')
plt.ylim(0,60)
plt.xlim(0,20)



ax2 = plt.subplot(gs[0, 5:9])

Xc,Yc=Aamp*100, Aeke
# X=Xc[Xc>=1];Y=Yc[Xc>=1];
# X=Xc[Xc>=np.nanmean(Xc)-np.nanstd(Xc)];Y=Yc[Xc>=np.nanmean(Xc)-np.nanstd(Xc)];
X=Xc[(Xc>=vamp)&(Wc>=vrad)];Y=Yc[(Xc>=vamp)&(Wc>=vrad)];


slope, intercept, r_value, p_value, std_err = stats.linregress(X,Y)
line = slope*(X)+intercept

plt.scatter(X,Y,c=colorx,edgecolors='white',linewidths=1,alpha=0.7)
plt.plot(X,line,color=colorx,alpha=0.5)

textstr = '\n'.join(('pearson r ={:.2f}'.format(r_value),
                    'p value ={:.2E}'.format(p_value)))


# plt.text(X.min(),Y.max()-0.5,textstr)
plt.title(textstr,fontweight='bold')
plt.grid()

plt.ylabel(r'EKE (cm$^2$ s$^{-2}$)',fontsize=fsz, fontweight='bold')
plt.xlabel('Amplitude (cm)',fontsize=fsz, fontweight='bold')
plt.ylim(0,Y.max())
plt.xlim(0,20)


ax3 = plt.subplot(gs[0, 10:14])

Xc,Yc=Aamp*100, Arad/1e3

# X=Xc[Xc>=1];Y=Yc[Xc>=1];
# X=Xc[Xc>=np.nanmean(Xc)-np.nanstd(Xc)];Y=Yc[Xc>=np.nanmean(Xc)-np.nanstd(Xc)];
X=Xc[(Xc>=vamp)&(Wc>=vrad)];Y=Yc[(Xc>=vamp)&(Wc>=vrad)];



slope, intercept, r_value, p_value, std_err = stats.linregress(X,Y)
line = slope*(X)+intercept

plt.scatter(X,Y,c=colorx,edgecolors='white',linewidths=1,alpha=0.7)
plt.plot(X,line,color=colorx,alpha=0.5)

textstr = '\n'.join(('pearson r ={:.2f}'.format(r_value),
                    'p value ={:.2E}'.format(p_value)))


# plt.text(X.min(),Y.max()-0.5,textstr)
plt.title(textstr,fontweight='bold')
plt.grid()

plt.ylabel(r'Radius (km)',fontsize=fsz, fontweight='bold')
plt.xlabel('Amplitude (cm)',fontsize=fsz, fontweight='bold')
plt.ylim(40,150)
plt.xlim(0,20)


ax4 = plt.subplot(gs[0,15:19])

xxc=Aamp*100
Xc,Yc=Arad/1e3,Apvel2*100;
# X=Xc[Xc>=1];Y=Yc[Xc>=1];
X=Xc[(xxc>=vamp)&(Wc>=vrad)];Y=Yc[(xxc>=vamp)&(Wc>=vrad)];

# X=Xc[Xc>=np.nanmean(Xc)-np.nanstd(Xc)];Y=Yc[Xc>=np.nanmean(Xc)-np.nanstd(Xc)];




slope, intercept, r_value, p_value, std_err = stats.linregress(X,Y)
line = slope*(X)+intercept

plt.scatter(X,Y,c=colorx,edgecolors='white',linewidths=1,alpha=0.7)
plt.plot(X,line,color=colorx,alpha=0.5)

textstr = '\n'.join(('pearson r ={:.2f}'.format(r_value),
                    'p value ={:.2E}'.format(p_value)))


# plt.text(X.min(),Y.max()-0.5,textstr)
plt.title(textstr,fontweight='bold')
plt.grid()

plt.ylabel(r'Propagating Speed (cm $s^{-1}$)',fontsize=fsz, fontweight='bold')
plt.xlabel('Radius (km)',fontsize=fsz, fontweight='bold')
plt.ylim(0,30)
plt.xlim(vrad,150)

plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/Scatter_AntiCyclone.png')
plt.close()


'''
Figura 2 - Histogramas
'''
Ra=Arad/1e3;
Rc=Crad/1e3;

plt.rcParams['font.weight']='bold'
plt.rcParams['font.size']=12
p=plt.figure(figsize=(16.57,  6.39))

gs = gridspec.GridSpec(6,11)
ax1 = plt.subplot(gs[0:3, 0:2])

Hc=-Camp*100
# cH=Hc.copy()[Hc<-1]
cH=Hc.copy()[(Hc<=-vamp)&(Rc>=40)]

# cH=Hc.copy()[Hc<np.nanmean(Hc)+np.nanstd(Hc)]


Ha=Aamp*100
# AA=Ah[Ha>1]

aH=Ha.copy()[(Ha>=vamp)&(Ra>=40)]
# aH=Ha.copy()[Ha>np.nanmean(Ha) -np.nanstd(Ha)]


plt.hist(cH,np.arange(-20,0,2),edgecolor='k',color='navy',alpha=0.7)
plt.hist(aH,np.arange(0,22,2),edgecolor='k',color='red',alpha=0.7)

plt.ylabel('Number of eddies',fontweight='bold')
plt.xlabel('Amplitude [cm]',fontweight='bold')

plt.xticks(np.linspace(-20,20,5),np.int_(np.abs(np.linspace(-20,20,5))))

ax2 = plt.subplot(gs[0:3, 3:5])

Ch=-Cvel*100
# CC=Ch[Hc<-1]
CC=Ch.copy()[(Hc<=-vamp)&(Rc>=40)]

# CC=Ch[Hc<np.nanmean(Hc)+np.nanstd(Hc)]



Ah=Avel*100
# AA=Ah[Ha>1]
AA=Ah.copy()[(Ha>=vamp)&(Ra>=40)]

# AA=Ah[Ha>np.nanmean(Ha) -np.nanstd(Ha)]

plt.hist(CC,np.arange(-40,0,4),edgecolor='k',color='navy',alpha=0.7)
plt.hist(AA,np.arange(0,44,4),edgecolor='k',color='red',alpha=0.7)
plt.ylabel('Number of eddies',fontweight='bold')
plt.xlabel(r'Swirling speed [cm s$^{-1}$]',fontweight='bold')

plt.xticks(np.linspace(-40,40,5),np.int_(np.abs(np.linspace(-40,40,5))))

ax3 = plt.subplot(gs[0:3, 6:8])

Ch=-Ceke
# CC=Ch[Hc<-1]
CC=Ch.copy()[(Hc<=-vamp)&(Rc>=40)]

# CC=Ch[Hc<np.nanmean(Hc)+np.nanstd(Hc)]


Ah=Aeke
# AA=Ah[Ha>1]
AA=Ah.copy()[(Ha>=vamp)&(Ra>=40)]

# AA=Ah[Ha>np.nanmean(Ha) -np.nanstd(Ha)]



plt.hist(CC,np.arange(-800,0,80),edgecolor='k',color='navy',alpha=0.7)
plt.hist(AA,np.arange(0,810,80),edgecolor='k',color='red',alpha=0.7)
plt.ylabel('Number of eddies',fontweight='bold')
plt.xlabel('EKE [cm$^{2}$ s$^{-2}$]',fontweight='bold')

plt.xticks(np.linspace(-800,800,5),np.int_(np.abs(np.linspace(-800,800,5))))


ax4 = plt.subplot(gs[0:3, 9:11])

Ch=-Crad/1e3
# CC=Ch[Hc<-1]
CC=Ch.copy()[(Hc<=-vamp)&(Rc>=40)]

# CC=Ch[Hc<np.nanmean(Hc)+np.nanstd(Hc)]


Ah=Arad/1e3
# AA=Ah[Ha>1]
AA=Ah.copy()[(Ha>=vamp)&(Ra>=40)]

# AA=Ah[Ha>np.nanmean(Ha) -np.nanstd(Ha)]

plt.hist(CC,np.arange(-200,0,20),edgecolor='k',color='navy',alpha=0.7)
plt.hist(AA,np.arange(0,220,20),edgecolor='k',color='red',alpha=0.7)
plt.ylabel('Number of eddies',fontweight='bold')
plt.xlabel('Radius [km]',fontweight='bold')

plt.xticks(np.linspace(-200,200,5),np.int_(np.abs(np.linspace(-200,200,5))))


ax5 = plt.subplot(gs[3:6, 0:2])


Ch=-Cpvel*100
# CC=Ch[Hc<-1]
CC=Ch.copy()[(Hc<=-vamp)&(Rc>=40)]

# CC=Ch[Hc<np.nanmean(Hc)+np.nanstd(Hc)]


Ah=Apvel*100
# AA=Ah[Ha>1]
AA=Ah.copy()[(Ha>=vamp)&(Ra>=40)]

# AA=Ah[Ha>np.nanmean(Ha) -np.nanstd(Ha)]


plt.hist(CC,np.arange(-50,0,4.5),edgecolor='k',color='navy',alpha=0.7)
plt.hist(AA,np.arange(0,50,4.5),edgecolor='k',color='red',alpha=0.7)
plt.ylabel('Number of eddies',fontweight='bold')
plt.xlabel(r'Propagating speed [cm s$^{-1}$]',fontweight='bold')

plt.xticks(np.linspace(-50,50,5),np.int_(np.abs(np.linspace(-50,50,5))))

ax6 = plt.subplot(gs[3:6, 3:5])



Ch=-Cdist
# CC=Ch[Hc<-1]
CC=Ch.copy()[(Hc<=-vamp)&(Rc>=40)]

# CC=Ch[Hc<np.nanmean(Hc)+np.nanstd(Hc)]


Ah=Adist
# AA=Ah[Ha>1]
AA=Ah.copy()[(Ha>=vamp)&(Ra>=40)]

# AA=Ah[Ha>np.nanmean(Ha) -np.nanstd(Ha)]


plt.hist(CC,np.arange(-2000,0,200),edgecolor='k',color='navy',alpha=0.7)
plt.hist(AA,np.arange(0,2200,200),edgecolor='k',color='red',alpha=0.7)
plt.ylabel('Number of eddies',fontweight='bold')
plt.xlabel('Distance [km]',fontweight='bold')


plt.xticks(np.linspace(-2000,2000,5),np.int_(np.abs(np.linspace(-2000,2000,5))))


ax7 = plt.subplot(gs[3:6, 6:8])


Ch=-Cdur
# CC=Ch[Hc<-1]
CC=Ch.copy()[(Hc<=-vamp)&(Rc>=40)]

# CC=Ch[Hc<np.nanmean(Hc)+np.nanstd(Hc)]


Ah=Adur
# AA=Ah[Ha>1]
AA=Ah.copy()[(Ha>=vamp)&(Ra>=40)]

# AA=Ah[Ha>np.nanmean(Ha) -np.nanstd(Ha)]


plt.hist(CC,np.arange(-250,0,25),edgecolor='k',color='navy',alpha=0.7)
plt.hist(AA,np.arange(0,275,25),edgecolor='k',color='red',alpha=0.7)
plt.ylabel('Number of eddies',fontweight='bold')
plt.xlabel('Duration [days]',fontweight='bold')

plt.xticks(np.linspace(-250,250,5),np.int_(np.abs(np.linspace(-250,250,5))))


plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/Total_histogram.png')
plt.close()









'Figuras 3a,3b,3c - Frequencies (Formation, Duration Monthly /Yearly)'
mostr=[u'Jan', u'Feb', u'Mar', u'Apr',u'May',u'Jun',u'Jul',u'Aug',u'Sep',u'Oct',u'Nov',u'Dec']


Bc=-Camp*100

cmo = pd.DatetimeIndex(Cdate).month
cfmo= pd.DatetimeIndex(Cfdate).month
cyea= pd.DatetimeIndex(Cdate).year

Cbox=np.zeros([12,len( np.arange(1993, Cdate[-1].year +1 ) )])


for mi,yi in zip(cmo,cyea):
    Cbox[mi-1,yi-1993]=len(cmo[(cmo==mi)&(cyea==yi)])

Cbox=pd.DataFrame(Cbox,columns=np.arange(1993, Cdate[-1].year +1 ),index=mostr)

Cbox.boxplot()


# cmo=cmo[Bc<np.nanmean(Bc) + np.nanstd(Bc)]
# cyea=cyea[Bc<np.nanmean(Bc) + np.nanstd(Bc)]

cyea=cyea[Bc<-vamp];cmo=cmo[Bc<-vamp];cfmo=cfmo[Bc<-vamp];

Ba=Aamp*100

amo = pd.DatetimeIndex(Adate).month
afmo= pd.DatetimeIndex(Afdate).month
ayea= pd.DatetimeIndex(Adate).year

amo=amo[Ba>vamp];ayea=ayea[Ba>vamp];afmo=afmo[Ba>vamp]

# amo=amo[Ba>np.nanmean(Ba)-np.nanstd(Ba)]
# ayea=ayea[Ba>np.nanmean(Ba)-np.nanstd(Ba)]



mrange=np.arange(0.5,13.5)
yrange=np.arange(cyea[0]-0.5,cyea[-1]+1 +0.5)

'''
CC
'''
p=plt.figure(figsize=(11.77,  3.86 ))

gs = gridspec.GridSpec(1,2)

ax1 = plt.subplot(gs[0, 0])


n=plt.hist(cmo,bins=mrange,color='navy',edgecolor='k',alpha=0.5)
# plt.plot(np.arange(1,13),n[0],'-',color='navy',alpha=1)

# sns.distplot(cmo,bins=12,color='navy',hist_kws={'edgecolor':'k'})
# sns.distplot(cfmo,bins=12,color='navy',hist_kws={'edgecolor':'k'})

ax1.set_xlabel('Months',fontweight='bold')
ax1.set_ylabel('Number of Eddies',fontweight='bold')

ax1.set_xlim(0,13)
ax1.set_ylim(n[0].min()-5,n[0].max()+5)
ax1.set_xticks(np.arange(1,13))
ax1.set_xticklabels(mostr)


ax2 = plt.subplot(gs[0,1])



n=plt.hist(cyea,bins=yrange,color='navy',edgecolor='k',alpha=0.5)
# plt.plot(np.arange(cyea[0],cyea[-1]+1),n[0],'-',color='navy',alpha=1)
# sns.distplot(cyea,bins=yrange,color='navy',hist_kws={'edgecolor':'k'})
ax2.set_xlabel('Years',fontweight='bold')

ax2.set_xlim(cyea[0]-1,cyea[-1]+1)
ax2.set_ylim(n[0].min()-5,n[0].max()+5)

ax2.set_xticks(np.arange(cyea[0],cyea[-1]+1)[::3])
ax2.set_xticklabels(np.arange(cyea[0],cyea[-1]+1)[::3])

plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/Total_frequency_Cyclones.png')
plt.close()


'''
AC
'''

p=plt.figure(figsize=(11.77,  3.86 ))
gs = gridspec.GridSpec(1,2)

ax1 = plt.subplot(gs[0, 0])

n=plt.hist(amo,bins=mrange,color='red',edgecolor='k',alpha=0.5)
# plt.plot(np.arange(1,13),n[0],'-',color='red',alpha=1)

# sns.distplot(amo,bins=12,color='red',hist_kws={'edgecolor':'k'})
# sns.distplot(afmo,bins=12,color='red',hist_kws={'edgecolor':'k'})

ax1.set_xlabel('Months',fontweight='bold')
ax1.set_ylabel('Number of Eddies',fontweight='bold')

ax1.set_xlim(0,13)
ax1.set_ylim(n[0].min()-5,n[0].max()+5)

ax1.set_xticks(np.arange(1,13))
ax1.set_xticklabels(mostr)


ax2 = plt.subplot(gs[0,1])




n=plt.hist(ayea,bins=yrange,color='red',edgecolor='k',alpha=0.5)
# plt.plot(np.arange(ayea[0],ayea[-1]+1),n[0],'-',color='red',alpha=1)

# sns.distplot(ayea,bins=yrange,color='red',hist_kws={'edgecolor':'k'})
ax2.set_xlabel('Years',fontweight='bold')

ax2.set_xlim(ayea[0]-1,ayea[-1]+1)
ax2.set_ylim(n[0].min()-5,n[0].max()+5)

ax2.set_xticks(np.arange(ayea[0],ayea[-1]+1)[::3])
ax2.set_xticklabels(np.arange(ayea[0],ayea[-1]+1)[::3])
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/Total_frequency_AntiCyclones.png')
plt.close()




'Figures 4 - Heatmaps of Dissipation, Formation, etc'



ltmn=-32;ltmx=-9;lnmn=-52;lnmx=-30.

resol=0.5
lonf=np.arange(lnmn,lnmx+resol,resol)
latf=np.arange(ltmn,ltmx+resol,resol)

# latf,lonf=np.meshgrid(latf,lonf)
lonf,latf=np.meshgrid(lonf,latf)

acmap=np.zeros(lonf.shape);
ccmap=np.zeros(lonf.shape);

'''
FORMATION REGIONS
'''
for ira in range(atrack.max()):
        arr=np.mean(anc['effective_radius'][:].data[atrack==ira])/1e3
        app=np.mean(anc['amplitude'][:].data[atrack==ira])*100

        if (arr>=40) & (app>=vamp):

            alx,aly=alon[atrack==ira][0],alat[atrack==ira][0]
            xind=np.argwhere(alx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(aly>=latf[:,0]).T[0][-1]

            acmap[yind,xind]+=1

for irc in range(ctrack.max()):
        crr=np.mean(cnc['effective_radius'][:].data[ctrack==irc])/1e3
        cpp=np.mean(cnc['amplitude'][:].data[ctrack==irc])*100

        if (crr>=40) & (cpp>=vamp):

            clx,cly=clon[ctrack==irc][0],clat[ctrack==irc][0]
            xind=np.argwhere(clx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(cly>=latf[:,0]).T[0][-1]
            ccmap[yind,xind]+=1



p=plt.figure(figsize=(8.21, 7.26))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=ltmn,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,2),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,2),labels=lbs,zorder=101,linewidth=0.1)

C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
# C.contourf(lonf,latf,((acmap+ ccmap)/np.nanmax(acmap + ccmap))*100,cmap='YlGnBu_r',latlon=True,levels=np.linspace(0,60,35))
plt.pcolormesh(lonf,latf,((acmap+ ccmap)/np.nanmax(acmap + ccmap))*100,cmap='hot_r',vmin=10,vmax=50);
cbar=plt.colorbar()
cbar.set_label('Formation [%]',fontweight='bold') 

C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)

plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/Total_formation.png')
plt.close()



'''
DISSIPATION REGIONS
'''
acdmap=np.zeros(lonf.shape);
ccdmap=np.zeros(lonf.shape);


for ira in range(atrack.max()):
        arr=np.mean(anc['effective_radius'][:].data[atrack==ira])/1e3
        app=np.mean(anc['amplitude'][:].data[atrack==ira])*100

        if (arr>=40) & (app>=vamp):

            alx,aly=alon[atrack==ira][-1],alat[atrack==ira][-1]
            xind=np.argwhere(alx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(aly>=latf[:,0]).T[0][-1]
            acdmap[yind,xind]+=1

for irc in range(ctrack.max()):
        crr=np.mean(cnc['effective_radius'][:].data[ctrack==irc])/1e3
        cpp=np.mean(cnc['amplitude'][:].data[ctrack==irc])*100

        if (crr>=40) & (cpp>=vamp):
            clx,cly=clon[ctrack==irc][-1],clat[ctrack==irc][-1]
            xind=np.argwhere(clx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(cly>=latf[:,0]).T[0][-1]
            ccdmap[yind,xind]+=1


p=plt.figure(figsize=(8.21, 7.26))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=ltmn,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,2),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,2),labels=lbs,zorder=101,linewidth=0.1)

C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
# C.contourf(lonf,latf,(acdmap + ccdmap)/np.nanmax(acdmap + ccdmap)*100,cmap='YlGnBu_r',latlon=True,levels=np.linspace(0,60,35))
plt.pcolormesh(lonf,latf,((acdmap+ ccdmap)/np.nanmax(acdmap + ccdmap))*100,cmap='hot_r',vmin=10,vmax=50);

cbar=plt.colorbar()
cbar.set_label('Dissipation [%]',fontweight='bold') 
C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/Total_dissipation.png')
plt.close()


'''
POLARITY BIRTH MAP
'''

'BOTH'
p=plt.figure(figsize=(8.21, 7.26))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=ltmn,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,2),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,2),labels=lbs,zorder=101,linewidth=0.1)

C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
C.contourf(lonf,latf,(acmap - ccmap)/np.nanmax(abs(acmap - ccmap))*100,cmap='seismic',latlon=True,levels=np.linspace(-100,100,35))
# plt.pcolormesh(lonf,latf,((acmap - ccmap)/np.nanmax(acmap - ccmap))*100,cmap='seismic',vmin=-100,vmax=100);

cbar=plt.colorbar(ticks=np.arange(-100,120,20))
cbar.set_label('Polarity Formation [%]',fontweight='bold')

C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/Total_polarity_birth.png')
plt.close()


'CYCLONES'
p=plt.figure(figsize=(8.21, 7.26))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=ltmn,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,2),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,2),labels=lbs,zorder=101,linewidth=0.1)

C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
# C.contourf(lonf,latf,(acmap - ccmap)/np.nanmax(abs(acmap - ccmap))*100,cmap='seismic',latlon=True,levels=np.linspace(-100,100,35))
plt.pcolormesh(lonf,latf,((ccmap)/np.nanmax(ccmap))*100,cmap='Blues',vmin=10,vmax=60);

cbar=plt.colorbar()
cbar.set_label(' Cyclone Formation [%]',fontweight='bold')

C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/cyclonic_birth.png')
plt.close()


'ANTICYCLONES'
p=plt.figure(figsize=(8.21, 7.26))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=ltmn,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,2),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,2),labels=lbs,zorder=101,linewidth=0.1)

C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
# C.contourf(lonf,latf,(acmap - ccmap)/np.nanmax(abs(acmap - ccmap))*100,cmap='seismic',latlon=True,levels=np.linspace(-100,100,35))
plt.pcolormesh(lonf,latf,((acmap)/np.nanmax(acmap))*100,cmap='Reds',vmin=10,vmax=60);

cbar=plt.colorbar()
cbar.set_label(' Anticyclone Formation [%]',fontweight='bold')

C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/anticyclonic_birth.png')
plt.close()


'''
RADIUS MAP
'''

radmap=np.zeros(lonf.shape);
adummy=np.zeros(lonf.shape);
cdummy=np.zeros(lonf.shape);


for ira in range(atrack.max()):
        arr=np.mean(anc['effective_radius'][:].data[atrack==ira])/1e3
        app=np.mean(anc['amplitude'][:].data[atrack==ira])*100

        if (arr>=40) & (app>=vamp):

            alx,aly=np.mean(alon[atrack==ira]),np.mean(alat[atrack==ira])

            xind=np.argwhere(alx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(aly>=latf[:,0]).T[0][-1]

            adummy[yind,xind]+=1
            radmap[yind,xind]+=Arad[ira]/1e3


for irc in range(ctrack.max()):
        crr=np.mean(cnc['effective_radius'][:].data[ctrack==irc])/1e3
        cpp=np.mean(cnc['amplitude'][:].data[ctrack==irc])*100

        if (crr>=40) & (cpp>=vamp):

            clx,cly=np.mean(clon[ctrack==irc]),np.mean(clat[ctrack==irc])
            xind=np.argwhere(clx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(cly>=latf[:,0]).T[0][-1]

            cdummy[yind,xind]+=1
            radmap[yind,xind]+=Crad[irc]/1e3



p=plt.figure(figsize=(9.47, 7.59))
C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=ltmn+1,urcrnrlat=ltmx-2,urcrnrlon=lnmx-1,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,2),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,2),labels=lbs,zorder=101,linewidth=0.1)

C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
# C.contourf(lonf,latf,(radmap)/(adummy+cdummy),cmap='rainbow',extend='both',latlon=True,levels=np.linspace(40,130,35))
# plt.pcolormesh(lonf,latf,(radmap)/(adummy+cdummy),cmap='rainbow',vmin=40,vmax=130);
plt.pcolormesh(lonf,latf,(radmap)/(adummy+cdummy),cmap=cm.cm.tempo,vmin=50,vmax=100);

cbar=plt.colorbar()
cbar.set_label('Radius [km]',fontweight='bold') 
C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)
plt.tight_layout()

plt.savefig('/home/igor/Documents/AntSimi/figures/Total_radius.png')
plt.close()



'''
AMP MAP
'''

ampmap=np.zeros(lonf.shape);
adummy=np.zeros(lonf.shape);
cdummy=np.zeros(lonf.shape);


for ira in range(atrack.max()):
        arr=np.mean(anc['effective_radius'][:].data[atrack==ira])/1e3
        app=np.mean(anc['amplitude'][:].data[atrack==ira])*100

        if (arr>=40) & (app>=vamp):

            alx,aly=np.mean(alon[atrack==ira]),np.mean(alat[atrack==ira])

            xind=np.argwhere(alx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(aly>=latf[:,0]).T[0][-1]

            adummy[yind,xind]+=1
            ampmap[yind,xind]+=Aamp[ira]*100


for irc in range(ctrack.max()):
        crr=np.mean(cnc['effective_radius'][:].data[ctrack==irc])/1e3
        cpp=np.mean(cnc['amplitude'][:].data[ctrack==irc])*100

        if (crr>=40) & (cpp>=vamp):

            clx,cly=np.mean(clon[ctrack==irc]),np.mean(clat[ctrack==irc])
            xind=np.argwhere(clx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(cly>=latf[:,0]).T[0][-1]

            cdummy[yind,xind]+=1
            ampmap[yind,xind]+=Camp[irc]*100





p=plt.figure(figsize=(9.47, 7.59))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=ltmn+1,urcrnrlat=ltmx-2,urcrnrlon=lnmx-1,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,2),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,2),labels=lbs,zorder=101,linewidth=0.1)

C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
# C.contourf(lonf,latf,(ampmap)/(adummy+cdummy),cmap='rainbow',extend='both',latlon=True,levels=np.linspace(vamp,10,35))
plt.pcolormesh(lonf,latf,(ampmap)/(adummy+cdummy),cmap=cm.cm.tempo,vmin=vamp+0.5,vmax=7);

cbar=plt.colorbar()
cbar.set_label('Amplitude [cm]',fontweight='bold') 
C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)
plt.tight_layout()

plt.savefig('/home/igor/Documents/AntSimi/figures/Total_amplitude.png')
plt.close()



'''
NON-LINEARITY
'''


nlmap=np.zeros(lonf.shape);
adummy=np.zeros(lonf.shape);
cdummy=np.zeros(lonf.shape);


for ira in range(atrack.max()):
        arr=np.mean(anc['effective_radius'][:].data[atrack==ira])/1e3
        app=np.mean(anc['amplitude'][:].data[atrack==ira])*100

        if (arr>=40) & (app>=vamp):

            alx,aly=np.mean(alon[atrack==ira]),np.mean(alat[atrack==ira])

            xind=np.argwhere(alx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(aly>=latf[:,0]).T[0][-1]

            adummy[yind,xind]+=1
            nlmap[yind,xind]+=(np.nanmax(anc['speed_average'][:][atrack==ira])/(Apvel2[ira]))


for irc in range(ctrack.max()):
        crr=np.mean(cnc['effective_radius'][:].data[ctrack==irc])/1e3
        cpp=np.mean(cnc['amplitude'][:].data[ctrack==irc])*100

        if (crr>=40) & (cpp>=vamp):

            clx,cly=np.mean(clon[ctrack==irc]),np.mean(clat[ctrack==irc])
            xind=np.argwhere(clx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(cly>=latf[:,0]).T[0][-1]

            cdummy[yind,xind]+=1
            nlmap[yind,xind]+=(np.nanmax(cnc['speed_average'][:][ctrack==irc])/(Cpvel2[irc]))




p=plt.figure(figsize=(9.47, 7.59))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=ltmn+1,urcrnrlat=ltmx-2,urcrnrlon=lnmx-1,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,2),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,2),labels=lbs,zorder=101,linewidth=0.1)

C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
# C.contourf(lonf,latf,(nlmap)/(adummy+cdummy),cmap='rainbow',extend='both',latlon=True,levels=np.linspace(0,6,35))
plt.pcolormesh(lonf,latf,(nlmap)/(adummy+cdummy),cmap=cm.cm.tempo,vmin=1,vmax=7);

cbar=plt.colorbar()
cbar.set_label('Advective Nonlinearity Parameter',fontweight='bold') 
C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/Total_nonlinear.png')
plt.close()
