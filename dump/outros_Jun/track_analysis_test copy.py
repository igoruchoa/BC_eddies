from netCDF4 import Dataset
from py_eddy_tracker.dataset.grid import RegularGridDataset, UnRegularGridDataset
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
from matplotlib.path import Path

'''
Plots efficiency and statiscal analysis of Anticyclones and cyclones
'''

'External Functions'
def smoo(self,nw):
    import scipy.signal as sg
    win = np.blackman(nw)
    win=win[None,:]*win[:,None]
    peso=sg.fftconvolve(np.ones(self.shape),win,mode='same')
    dataout=sg.fftconvolve(self,win,mode='same')/peso

    return dataout




def geocut(xb,yb,var,ltmn,lnmn,ltmx,lnmx):
    if len(xb.shape)>1:
        indx=np.argwhere((xb[0,:]>=lnmn)&(xb[0,:]<=lnmx)).T[0]
        indy=np.argwhere((yb[:,0]>=ltmn)&(yb[:,0]<=ltmx)).T[0]
        xnew=xb[indy,:][:,indx]
        ynew=yb[indy,:][:,indx]

    else:
        indx=np.argwhere((xb>=lnmn)&(xb<=lnmx)).T[0]
        indy=np.argwhere((yb>=ltmn)&(yb<=ltmx)).T[0]
        xnew=xb[indx]
        ynew=yb[indy]

    varnew=var[indy,:][:,indx]
    return xnew,ynew,varnew


# 'Loading Data'
# var=Dataset('/home/igor/Documents/AntSimi/altdata/clim/sla_climatology.nc')

# lonc,latc=var['longitude'][:]-360,var['latitude'][:];
# slac=var['sla'][:]

dt0=datetime(1950,1,1)
# dt=[dt0 + timedelta(days=np.int(dia)) for dia in var['time']]

'Loading Tracks'
anc=Dataset('/home/igor/Documents/META_eddy_detection/META3.1exp_DT_allsat_Anticyclonic_long_19930101_20200307.nc')
cnc=Dataset('/home/igor/Documents/META_eddy_detection/META3.1exp_DT_allsat_Cyclonic_long_19930101_20200307.nc')


# anc=Dataset('/home/igor/Documents/AntSimi/tracking/Anticyclonic.nc')
# cnc=Dataset('/home/igor/Documents/AntSimi/tracking/Cyclonic.nc')



'Delimited Area'

import pandas as pd
# dd=pd.read_csv('/home/igor/Documents/AntSimi/data/buffer_200m_6rad.csv')
dd=pd.read_csv('/home/igor/Documents/AntSimi/data/bat/buffer_200m_200km.csv')
ddlon,ddlat=np.array(dd['lon']),np.array(dd['lat'])
ddseg=Path(np.array([ddlon,ddlat]).T)


# ctrack=pd.Series(ctrack.data)
# ctrack=ctrack.data;

'Topography'
ltmn=-29;ltmx=-14;lnmn=-52;lnmx=-35.

topo=Dataset('/home/igor/Documents/Extras_Works/Sumida/HYCOM/bat/etopo1_bedrock.nc')
xb=topo['lon'][:];yb=topo['lat'][:];zb=topo['Band1'][:];
zb=smoo(zb,20)
xb,yb,zb=geocut(xb,yb,zb,ltmn,lnmn,ltmx,lnmx);
xxb,yyb=np.meshgrid(xb,yb);


'Bathymatric Path'
c=plt.contourf(xb,yb,zb,levels=[-200,np.inf],colors='k')
plt.close()
segs=[]

for iss in range(len(c.allsegs)):
    # col=c.collections[iss].get_paths()
    col=c.allsegs[iss]

    lens=list(map(len,col))
    indl=np.argmax(lens)
    segs.append(col[indl])

segs=np.vstack(segs);
pseg=Path(segs)


atrack=anc['track'][:]
aobs=anc['observation_number'][:]
alon,alat=anc['longitude'][:],anc['latitude'][:]
atrack=np.int_(atrack.data)
alon=np.where(alon>180,alon-360,alon)
# atrack=pd.Series(atrack.data)
# atrack=atrack.data;

ctrack=cnc['track'][:]
cobs=cnc['observation_number'][:]
clon,clat=cnc['longitude'][:],cnc['latitude'][:]
ctrack=np.int_(ctrack.data)
clon=np.where(clon>180,clon-360,clon)

ltmn=-29;ltmx=-14;lnmn=-52;lnmx=-35.


'Condition 1 - Excluding NH++'
condc=(clon>lnmn-10)&(clon<lnmx+10)&(clat>ltmn-10)&(clat<ltmx+10) #esse é o certo
conda=(alon>lnmn-10)&(alon<lnmx+10)&(alat>ltmn-10)&(alat<ltmx+10)



cchunck=ctrack[condc]
achunck=atrack[conda]
ctrack[~condc]=-1
atrack[~conda]=-1


'Condition 1b - Filtering NH++'
xxc=(ctrack>=0)&(cobs==0)
condc1=(ctrack[condc]>=0)&(cobs[condc]==0)&(clon[condc]>lnmn)&(clon[condc]<lnmx)&(clat[condc]>ltmn)&(clat[condc]<ltmx)
ics=ctrack[condc][~condc1] #indices iniciais que precisam ser retirados com  o resto de suas tracks

# condc=(clon[xxc]>lnmn)&(clon[xxc]<lnmx)&(clat[xxc]>ltmn)&(clat[xxc]<ltmx)
# ics=ctrack[xxc][~condc] #indices iniciais que precisam ser retirados com  o resto de suas tracks


def Completer(track,ind):
    xtrack=track.copy()
    iii=0
    for iy in ind:

        print('{0}%'.format((iii/len(ind))*100),end='\r')
        xtrack=np.where(xtrack==iy,-1,xtrack)
        iii+=1
    return xtrack    

ctrack1=Completer(ctrack,ics)


xxa=(atrack>=0)&(aobs==0)
conda=(alon[xxa]>lnmn)&(alon[xxa]<lnmx)&(alat[xxa]>ltmn)&(alat[xxa]<ltmx)
ias=atrack[xxa][~conda] #indices que precisam ser retirados com  o resto de suas tracks

atrack1=Completer(atrack,ias)



'Condition 2 - Bathtymetry'
xxc=(ctrack1>=0)&(cobs==0)
cbol=[pseg.contains_point((clon[iy],clat[iy])) for iy in np.where(xxc)[0] ]
ics=ctrack1[xxc][np.array(cbol)]

ctrack2=Completer(ctrack1,ics)


xxa=(atrack1>=0)&(aobs==0)
abol=[pseg.contains_point((alon[iy],alat[iy])) for iy in np.where(xxa)[0] ]
ias=atrack1[xxa][np.array(abol)]

atrack2=Completer(atrack1,ias)


'Condition 3 - Buffer zone'
xxc=(ctrack2>=0)&(cobs==0)
cuni=[ddseg.contains_point((clon[iy],clat[iy])) for iy in np.where(xxc)[0] ]
ics=ctrack2[xxc][~np.array(cuni)]

ctrack3=Completer(ctrack2,ics)


xxa=(atrack2>=0)&(aobs==0)
auni=[ddseg.contains_point((alon[iy],alat[iy])) for iy in np.where(xxa)[0] ]
ias=atrack2[xxa][~np.array(auni)]

atrack3=Completer(atrack2,ias)


pd.Series.to_pickle(atrack3,'/home/igor/Documents/META_eddy_detection/atrack_final.pkl')
pd.Series.to_pickle(ctrack3,'/home/igor/Documents/META_eddy_detection/ctrack_final.pkl')

atrack=atrack3.copy()
ctrack=ctrack3.copy()
del atrack3,atrack1,atrack2,ctrack3,ctrack1,ctrack2


# atrack3=pd.read_pickle('/home/igor/Documents/META_eddy_detection/atrack_final.pkl')
# ctrack3=pd.read_pickle('/home/igor/Documents/META_eddy_detection/ctrack_final.pkl')

# atrack=atrack3.copy()
# ctrack=ctrack3.copy()
# del atrack3,ctrack3
'################################################################################'



'Filtering final'
# ctrack=pd.Series()
# atrack=np.float_(atrack);

# ctrack[ctrack<0]=np.nan;
# atrack[atrack<0]=np.nan;
aunik=np.unique(atrack[atrack>=0])
cunik=np.unique(ctrack[ctrack>=0])





'Snapshot Tracked'
ref=anc['time'][:][0]
at=np.int_((anc['time'][:]-ref))
ct=np.int_((cnc['time'][:]-ref))


'Anticyclonic parameters'
Atime=anc['time'][:].data
Aamp=pd.Series(anc['amplitude'][:]);Aamp=Aamp.groupby(atrack).mean();
Avel=pd.Series(anc['speed_average'][:]);Avel=Avel.groupby(atrack).mean();
Arad=pd.Series(anc['effective_radius'][:]);Arad=Arad.groupby(atrack).mean();
# AyContour=pd.DataFrame(anc['effective_contour_latitude'][:]);
# AxContour=pd.DataFrame(anc['effective_contour_longitude'][:])-360;

Aeke=0.5*((Avel*100)**2)

Aamp=Aamp[1:]
Avel=Avel[1:]
Arad=Arad[1:]
Aeke=Aeke[1:]


Adist=[];Adur=[];Apvel=[];Apvel2=[];Apvel3=[];
Adate=[];Afdate=[];

iii=0
for ai in aunik:
# for ai in np.unique(atrack):
# for ai in range(atrack.max()+1):
    print('{0}%'.format((iii/len(aunik))*100),end='\r')
    iii+=1

    #duration
    dur=len(atrack[atrack==ai])
    #distance
    distx=sw.dist(lon=alon[atrack==ai],lat=alat[atrack==ai])[0]
    disty=sw.dist(lon=[alon[atrack==ai][0],alon[atrack==ai][-1]],
                    lat=[alat[atrack==ai][0],alat[atrack==ai][-1]])[0]
    # distz=sw.dist(lon=[AxContour.iloc[atrack==ai,:].mean(1).iloc[0],AxContour.iloc[atrack==ai,:].mean(1).iloc[-1]],
    #                 lat=[AyContour.iloc[atrack==ai,:].mean(1).iloc[0],AyContour.iloc[atrack==ai,:].mean(1).iloc[-1]])[0]


    k=np.insert(distx,0,0)
    dist=np.cumsum(k)
    #propagation speed
    pvel=(dist[-1]*1e3)/(dur*86400)
    pvel2=(disty[0]*1e3/(dur*86400))
    # pvel3=(distz[0]*1e3/(dur*86400))
    #time
    time=Atime[atrack==ai];

    date=dt0+ timedelta(days=np.int(time.mean()));
    fdate=dt0 + timedelta(days=np.int(time[0]));

    #Concatenating
    
    # Adist.append(dist[-1]);
    Adist.append(disty[0]);

    Adur.append(dur);
    Apvel.append(pvel*100);
    Apvel2.append(pvel2*100);
    # Apvel3.append(pvel3*100);
    Adate.append(date);
    Afdate.append(fdate);

Adist=np.array(Adist);
Adur=np.array(Adur);
Apvel=np.array(Apvel);
Apvel2=np.array(Apvel2);
# Apvel3=np.array(Apvel3);

'Cyclonic parameters'
Ctime=cnc['time'][:].data
Camp=pd.Series(cnc['amplitude'][:]);Camp=Camp.groupby(ctrack).mean();
Cvel=pd.Series(cnc['speed_average'][:]);Cvel=Cvel.groupby(ctrack).mean();
Crad=pd.Series(cnc['effective_radius'][:]);Crad=Crad.groupby(ctrack).mean();
# CyContour=pd.DataFrame(cnc['effective_contour_latitude'][:]);
# CxContour=pd.DataFrame(cnc['effective_contour_longitude'][:])-360;


Ceke=0.5*((Cvel*100)**2)

Camp=Camp[1:]
Cvel=Cvel[1:]
Crad=Crad[1:]
Ceke=Ceke[1:]

Cdist=[];Cdur=[];Cpvel=[];Cpvel2=[];Cpvel3=[];
Cdate=[];Cfdate=[];


iii=0
for ci in cunik:
# for ci in range(ctrack.max()+1):
# for ci in np.unique(ctrack):

    print('{0}%'.format((iii/len(cunik))*100),end='\r')
    iii+=1
    dur=len(ctrack[ctrack==ci])

    distx=sw.dist(lon=clon[ctrack==ci],lat=clat[ctrack==ci])[0]
    disty=sw.dist(lon=[clon[ctrack==ci][0],clon[ctrack==ci][-1]],
                    lat=[clat[ctrack==ci][0],clat[ctrack==ci][-1]])[0]
    # distz=sw.dist(lon=[CxContour.iloc[ctrack==ci,:].mean(1).iloc[0],CxContour.iloc[ctrack==ci,:].mean(1).iloc[-1]],
                    # lat=[CyContour.iloc[ctrack==ci,:].mean(1).iloc[0],CyContour.iloc[ctrack==ci,:].mean(1).iloc[-1]])[0]




    k=np.insert(distx,0,0)
    dist=np.cumsum(k)

    pvel=(dist[-1]*1e3)/(dur*86400)
    pvel2=(disty[0]*1e3/(dur*86400))
    # pvel3=(distz[0]*1e3/(dur*86400))

    time=Ctime[ctrack==ci];

    date=dt0+ timedelta(days=np.int(time.mean()));
    fdate=dt0 + timedelta(days=np.int(time[0]));

    # Cdist.append(dist[-1]);
    Cdist.append(disty)
    Cdur.append(dur);
    Cpvel.append(pvel*100);
    Cpvel2.append(pvel2*100);
    # Cpvel3.append(pvel3*100);
    Cdate.append(date);
    Cfdate.append(fdate);


Cdist=np.array(Cdist)
Cdur=np.array(Cdur)
Cpvel=np.array(Cpvel)
Cpvel2=np.array(Cpvel2)
# Cpvel3=np.array(Cpvel3)


#Contratining/Filtering by amplitude and radius
vamp=1 #cm
vrad=40 #km

#masks 
# ccond=(Camp>=vamp/100)&(Crad>=vrad*1e3) # m2cm / m2km
# acond=(Aamp>=vamp/100)&(Arad>=vrad*1e3)

ccond=(Camp==Camp) # m2cm / m2km
acond=(Aamp==Aamp)

'''
#plot the daily images with detections
# for ip in range(slac.shape[0]):
for ip in range(365):
    p=plt.figure(figsize=(4.5 , 4.5))
    plt.contourf(lonc,latc,slac[ip],cmap='seismic',levels=np.linspace(-0.2,0.2,50))
    plt.contour(lonc,latc,slac[ip],colors='dimgray',linewidths=0.5,levels=np.linspace(-0.2,0.2,30))


    # asx=np.unique(atrack[at==ip])
    # for yi in asx[asx>0]:
    #     if yi in aunik:
    #         plt.plot(alon[atrack==yi],alat[atrack==yi],color='IndianRed',marker='o',markersize=3)
    #         plt.text(alon[atrack==yi],alat[atrack==yi],np.int(yi),color='IndianRed',fontweight='bold')


    # csx=np.unique(ctrack[ct==ip])
    # for yi in csx[csx>0]:
    #     if yi in cunik:
    #         plt.plot(clon[ctrack==yi][0],clat[ctrack==yi][0],color='navy',marker='o',markersize=3)
    #         plt.text(clon[ctrack==yi][0],clat[ctrack==yi][0],np.int(yi),color='navy',fontweight='bold')

    [plt.plot(xll,yll,'ro',markersize=3) for xll,yll in zip(alon[at==ip][atrack[at==ip]>0],alat[at==ip][atrack[at==ip]>0]) ]
    [plt.text(xll,yll,np.int(tll),color='r',fontweight='bold') for xll,yll,tll in zip(alon[at==ip][atrack[at==ip]>0],alat[at==ip][atrack[at==ip]>0],atrack[at==ip][atrack[at==ip]>0]) ]

    [plt.plot(xll,yll,'bo',markersize=3) for xll,yll in zip(clon[ct==ip][ctrack[ct==ip]>0],clat[ct==ip][ctrack[ct==ip]>0]) ]
    [plt.text(xll,yll,np.int(tll),color='b',fontweight='bold') for xll,yll,tll in zip(clon[ct==ip][ctrack[ct==ip]>0],clat[ct==ip][ctrack[ct==ip]>0],ctrack[ct==ip][ctrack[ct==ip]>0]) ]
    
    plt.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-200,np.inf])
    plt.contourf(xxb,yyb,zb,colors='green',alpha=1,levels=[0,np.inf])

    plt.ylim(-22,-14)
    plt.xlim(-41,-33)

    # plt.tight_layout()
    plt.savefig('/home/igor/Documents/META_eddy_detection/maps/{0:04d}_track.png'.format(ip+1))
    plt.close()
'''


'Plotting Eddy Formation'
plt.rcParams['font.size']=20
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.weight']='bold'

p=plt.figure(figsize=(9.26, 7.96))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmn+2,llcrnrlat=ltmn,urcrnrlat=ltmx,urcrnrlon=lnmx+0.2,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,4),labels=lbs,zorder=101,linewidth=0.1)

plt.scatter(alon[atrack>0][aobs[atrack>0]==0],alat[atrack>0][aobs[atrack>0]==0],alpha=0.5,s=4,color='r')
# [C.plot(alon[atrack==iy],alat[atrack==iy],'r',linewidth=0.5) for iy in range(atrack.max())]
# [C.plot(alon[atrack==iy][0],alat[atrack==iy][0],'r.',markersize=5,alpha=0.5) for iy in aunik[acond]]

plt.scatter(clon[ctrack>0][cobs[ctrack>0]==0],clat[ctrack>0][cobs[ctrack>0]==0],alpha=0.5,s=4,color='b')
# [C.plot(clon[ctrack==iy],clat[ctrack==iy],'b',linewidth=0.5) for iy in range(ctrack.max())]
# [C.plot(clon[ctrack==iy][0],clat[ctrack==iy][0],'b.',markersize=5,alpha=0.5) for iy in cunik[ccond]]

plt.contourf(xb,yb,zb,levels=[-200,np.inf],colors='dimgray')

cq=plt.contour(xb,yb,zb,levels=[-200],colors='k')
# plt.clabel(cq,fmt='%.d')

plt.tight_layout()
plt.savefig('/home/igor/Dropbox/META_eddydetection_images/formation_points.png')
plt.close()



'Plotting Eddy Dissipation'

p=plt.figure(figsize=(9.26, 7.96))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmn+2,llcrnrlat=ltmn,urcrnrlat=ltmx,urcrnrlon=lnmx+0.2,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1,zorder=110);
C.fillcontinents(color='#fecf78',zorder=109)
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,zorder=111,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,4),labels=lbs,zorder=112,linewidth=0.1)

# [C.plot(alon[atrack==iy],alat[atrack==iy],'r',linewidth=0.5) for iy in range(atrack.max())]
[C.plot(alon[atrack==iy][-1],alat[atrack==iy][-1],'r.',markersize=5,alpha=0.5) for iy in aunik[acond]]

# [C.plot(clon[ctrack==iy],clat[ctrack==iy],'b',linewidth=0.5) for iy in range(ctrack.max())]
[C.plot(clon[ctrack==iy][-1],clat[ctrack==iy][-1],'b.',markersize=5,alpha=0.5) for iy in cunik[ccond]]

plt.contourf(xb,yb,zb,levels=[-200,np.inf],colors='dimgray',zorder=101)

cq=plt.contour(xb,yb,zb,levels=[-200],colors='k',zorder=105)

# xc=plt.clabel(cq,fmt='%.d')


plt.tight_layout()
plt.savefig('/home/igor/Dropbox/META_eddydetection_images/dissipation_points.png')
plt.close()



#Saving table of statistics
'Table'
ctt=pd.Series({'amp':np.nanmean(Camp[ccond]),'amp_min':np.nanmin(Camp[ccond]),'amp_max':np.nanmax(Camp[ccond]),'amp_std':np.nanstd(Camp[ccond]),
                'vel':np.nanmean(Cvel[ccond]),'vel_min':np.nanmin(Cvel[ccond]),'vel_max':np.nanmax(Cvel[ccond]),'vel_std':np.nanstd(Cvel[ccond]),
                'rad':np.nanmean(Crad[ccond]),'rad_min':np.nanmin(Crad[ccond]),'rad_max':np.nanmax(Crad[ccond]),'rad_std':np.nanstd(Crad[ccond]),
                'eke':np.nanmean(Ceke[ccond]),'eke_min':np.nanmin(Ceke[ccond]),'eke_max':np.nanmax(Ceke[ccond]),'eke_std':np.nanstd(Ceke[ccond]),
                'dist':np.nanmean(Cdist[ccond]),'dist_min':np.nanmin(Cdist[ccond]),'dist_max':np.nanmax(Cdist[ccond]),'dist_std':np.nanstd(Cdist[ccond]),
                'dur':np.nanmean(Cdur[ccond]),'dur_min':np.nanmin(Cdur[ccond]),'dur_max':np.nanmax(Cdur[ccond]),'dur_std':np.nanstd(Cdur[ccond]),
                'pvel':np.nanmean(Cpvel[ccond]),'pvel_min':np.nanmin(Cpvel[ccond]),'pvel_max':np.nanmax(Cpvel[ccond]),'pvel_std':np.nanstd(Cpvel[ccond]),
                'pvel2':np.nanmean(Cpvel2[ccond]),'pvel2_min':np.nanmin(Cpvel2[ccond]),'pvel2_max':np.nanmax(Cpvel2[ccond]),'pvel2_std':np.nanstd(Cpvel2[ccond])})
                # 'pvel3':np.nanmean(Cpvel3[ccond]),'pvel2_min':np.nanmin(Cpvel3[ccond]),'pvel3_max':np.nanmax(Cpvel3[ccond]),'pvel3_std':np.nanstd(Cpvel3[ccond])})


att=pd.Series({'amp':np.nanmean(Aamp[acond]),'amp_min':np.nanmin(Aamp[acond]),'amp_max':np.nanmax(Aamp[acond]),'amp_std':np.nanstd(Aamp[acond]),
                'vel':np.nanmean(Avel[acond]),'vel_min':np.nanmin(Avel[acond]),'vel_max':np.nanmax(Avel[acond]),'vel_std':np.nanstd(Avel[acond]),
                'rad':np.nanmean(Arad[acond]),'rad_min':np.nanmin(Arad[acond]),'rad_max':np.nanmax(Arad[acond]),'rad_std':np.nanstd(Arad[acond]),
                'eke':np.nanmean(Aeke[acond]),'eke_min':np.nanmin(Aeke[acond]),'eke_max':np.nanmax(Aeke[acond]),'eke_std':np.nanstd(Aeke[acond]),
                'dist':np.nanmean(Adist[acond]),'dist_min':np.nanmin(Adist[acond]),'dist_max':np.nanmax(Adist[acond]),'dist_std':np.nanstd(Adist[acond]),
                'dur':np.nanmean(Adur[acond]),'dur_min':np.nanmin(Adur[acond]),'dur_max':np.nanmax(Adur[acond]),'dur_std':np.nanstd(Adur[acond]),
                'pvel':np.nanmean(Apvel[acond]),'pvel_min':np.nanmin(Apvel[acond]),'pvel_max':np.nanmax(Apvel[acond]),'pvel_std':np.nanstd(Apvel[acond]),
                'pvel2':np.nanmean(Apvel2[acond]),'pvel2_min':np.nanmin(Apvel2[acond]),'pvel2_max':np.nanmax(Apvel2[acond]),'pvel2_std':np.nanstd(Apvel2[acond])})
                # 'pvel3':np.nanmean(Apvel3[acond]),'pvel3_min':np.nanmin(Apvel3[acond]),'pvel3_max':np.nanmax(Apvel3[acond]),'pvel3_std':np.nanstd(Apvel3[acond])})



savep='/home/igor/Dropbox/META_eddydetection_images/'
ctt.to_csv(savep+ 'cyclone_prop.dat',sep='\t',mode='w',float_format='%.6f',index=True,header=True,na_rep='NaN')
att.to_csv(savep + 'anticyclone_prop.dat',sep='\t',mode='w',float_format='%.6f',index=True,header=True,na_rep='NaN')



#Number of eddies
print('Number of Eddies: {0}'.format(aunik[acond].shape[0] + cunik[ccond].shape[0]))



print('Number of AC: {0}'.format(aunik[acond].shape[0]))
anumber=aunik[acond].shape[0]

print('Number of CC: {0}'.format(cunik[ccond].shape[0]))
cnumber=cunik[ccond].shape[0]




#Plots
'''
Figure 1a & 1b - Regression Plots
CC
'''

Wc=Crad/1e3

colorx='navy'
fsz=20;

pc=plt.figure(figsize=(18.53,  5.46))

gs = gridspec.GridSpec(1,14)
ax1 = plt.subplot(gs[0, 0:4])

Xc,Yc=-Camp*100,Cvel*100
# X=Xc[Xc<=-1];Y=Yc[Xc<=-1];
X=Xc[(Xc<=-vamp)&(Wc>=vrad)];Y=Yc[(Xc<=-vamp)&(Wc>=vrad)];


slope, intercept, r_value, p_value, std_err = stats.linregress(X,Y)
line = slope*(X)+intercept

plt.scatter(X,Y,c=colorx,edgecolors='white',linewidths=1,alpha=0.7)
plt.plot(X,line,color=colorx,alpha=0.5)

# textstr = '\n'.join(('pearson r ={:.2f}'.format(r_value),
#                     'p value ={0}'.format(p_value)))
                    
textstr = ('pearson r ={:.2f}'.format(r_value))

plt.title(textstr,fontweight='bold',fontsize=fsz-2)

plt.grid()

plt.ylabel(r'Swirling Speed (cm s$^{-1}$)',fontsize=fsz, fontweight='bold')
plt.xlabel('Amplitude (cm)',fontsize=fsz, fontweight='bold')
plt.ylim(0,100)
plt.xlim(-60,0)



ax2 = plt.subplot(gs[0, 5:9])

Xc,Yc=-Camp*100,Ceke
# X=Xc[Xc<=-1];Y=Yc[Xc<=-1];
X=Xc[(Xc<=-vamp)&(Wc>=vrad)];Y=Yc[(Xc<=-vamp)&(Wc>=vrad)];

slope, intercept, r_value, p_value, std_err = stats.linregress(X,Y)
line = slope*(X)+intercept

plt.scatter(X,Y,c=colorx,edgecolors='white',linewidths=1,alpha=0.7)
plt.plot(X,line,color=colorx,alpha=0.5)

# textstr = '\n'.join(('pearson r ={:.2f}'.format(r_value),
#                     'p value ={0}'.format(p_value)))
                    
textstr = ('pearson r ={:.2f}'.format(r_value))

plt.title(textstr,fontweight='bold',fontsize=fsz-2)
plt.grid()

plt.ylabel(r'EKE (cm$^2$ s$^{-2}$)',fontsize=fsz, fontweight='bold')
plt.xlabel('Amplitude (cm)',fontsize=fsz, fontweight='bold')
plt.ylim(0,3000)
plt.xlim(-60,0)


ax3 = plt.subplot(gs[0, 10:14])

Xc,Yc=-Camp*100, Crad/1e3
# X=Xc[Xc<=-1];Y=Yc[Xc<=-1];
# X=Xc[(Xc<=-1 )&( Yc>=40)];Y=Yc[(Xc<=-1 )&(Yc>=40)];
X=Xc[(Xc<=-vamp)&(Wc>=vrad)];Y=Yc[(Xc<=-vamp)&(Wc>=vrad)];




slope, intercept, r_value, p_value, std_err = stats.linregress(X,Y)
line = slope*(X)+intercept

plt.scatter(X,Y,c=colorx,edgecolors='white',linewidths=1,alpha=0.7)
plt.plot(X,line,color=colorx,alpha=0.5)

# textstr = '\n'.join(('pearson r ={:.2f}'.format(r_value),
#                     'p value ={0}'.format(p_value)))
                    
textstr = ('pearson r ={:.2f}'.format(r_value))


# plt.text(X.min(),Y.max()-0.5,textstr)
plt.title(textstr,fontweight='bold',fontsize=fsz-2)
plt.grid()

plt.ylabel(r'Radius (km)',fontsize=fsz, fontweight='bold')
plt.xlabel('Amplitude (cm)',fontsize=fsz, fontweight='bold')
plt.ylim(30,200)
plt.xlim(-60,0)


# ax4 = plt.subplot(gs[0,15:19])

# xxc=-Camp*100
# Xc,Yc=-Crad/1e3, Cpvel3
# # X=Xc[Xc<=-1];Y=Yc[Xc<=-1];
# X=Xc[(xxc<=-vamp)&(Wc>=vrad)];Y=Yc[(xxc<=-vamp)&(Wc>=vrad)];



# slope, intercept, r_value, p_value, std_err = stats.linregress(X,Y)
# line = slope*(X)+intercept

# plt.scatter(X,Y,c=colorx,edgecolors='white',linewidths=1,alpha=0.7)
# plt.plot(X,line,color=colorx,alpha=0.5)

# # textstr = '\n'.join(('pearson r ={:.2f}'.format(r_value),
# #                     'p value ={0}'.format(p_value)))
                    
# textstr = ('pearson r ={:.2f}'.format(r_value))


# # plt.text(X.min(),Y.max()-0.5,textstr)
# plt.title(textstr,fontweight='bold',fontsize=fsz-2)
# plt.grid()

# plt.ylabel(r'Propagating Speed (cm $s^{-1}$)',fontsize=fsz, fontweight='bold')
# plt.xlabel('Radius (km)',fontsize=fsz, fontweight='bold')
# plt.ylim(0,30)
# plt.xlim(-150,-30)

pc.subplots_adjust(top=0.9,
bottom=0.17,
left=0.06,
right=0.97,
hspace=0.21,
wspace=0.55)
plt.tight_layout()

plt.savefig('/home/igor/Dropbox/META_eddydetection_images/Scatter_Cyclone.png')
plt.close()



'''
AC
'''

Wc=Arad/1e3
colorx='red'
fsz=20;

pc=plt.figure(figsize=(18.53,  5.46))



gs = gridspec.GridSpec(1,14)
ax1 = plt.subplot(gs[0, 0:4])

Xc,Yc=Aamp*100,Avel*100

# X=Xc[Xc>=1];Y=Yc[Xc>=1];
# X=Xc[Xc>=np.nanmean(Xc)-np.nanstd(Xc)];Y=Yc[Xc>=np.nanmean(Xc)-np.nanstd(Xc)];
X=Xc[(Xc>=vamp)&(Wc>=vrad)];Y=Yc[(Xc>=vamp)&(Wc>=vrad)];

slope, intercept, r_value, p_value, std_err = stats.linregress(X,Y)
line = slope*(X)+intercept

plt.scatter(X,Y,c=colorx,edgecolors='white',linewidths=1,alpha=0.7)
plt.plot(X,line,color=colorx,alpha=0.5)

# textstr = '\n'.join(('pearson r ={:.2f}'.format(r_value),
#                     'p value ={0}'.format(p_value)))
                    
textstr = ('pearson r ={:.2f}'.format(r_value))


# plt.text(X.min(),Y.max()-0.5,textstr)
plt.title(textstr,fontweight='bold',fontsize=fsz-2)
plt.grid()

plt.ylabel(r'Swirling Speed (cm s$^{-1}$)',fontsize=fsz, fontweight='bold')
plt.xlabel('Amplitude (cm)',fontsize=fsz, fontweight='bold')
plt.ylim(0,100)
plt.xlim(0,60)



ax2 = plt.subplot(gs[0, 5:9])

Xc,Yc=Aamp*100, Aeke
# X=Xc[Xc>=1];Y=Yc[Xc>=1];
# X=Xc[Xc>=np.nanmean(Xc)-np.nanstd(Xc)];Y=Yc[Xc>=np.nanmean(Xc)-np.nanstd(Xc)];
X=Xc[(Xc>=vamp)&(Wc>=vrad)];Y=Yc[(Xc>=vamp)&(Wc>=vrad)];


slope, intercept, r_value, p_value, std_err = stats.linregress(X,Y)
line = slope*(X)+intercept

plt.scatter(X,Y,c=colorx,edgecolors='white',linewidths=1,alpha=0.7)
plt.plot(X,line,color=colorx,alpha=0.5)

# textstr = '\n'.join(('pearson r ={:.2f}'.format(r_value),
#                     'p value ={0}'.format(p_value)))
                    
textstr = ('pearson r ={:.2f}'.format(r_value))


# plt.text(X.min(),Y.max()-0.5,textstr)
plt.title(textstr,fontweight='bold',fontsize=fsz-2)
plt.grid()

plt.ylabel(r'EKE (cm$^2$ s$^{-2}$)',fontsize=fsz, fontweight='bold')
plt.xlabel('Amplitude (cm)',fontsize=fsz, fontweight='bold')
plt.ylim(0,3000)
plt.xlim(0,60)


ax3 = plt.subplot(gs[0, 10:14])

Xc,Yc=Aamp*100, Arad/1e3

# X=Xc[Xc>=1];Y=Yc[Xc>=1];
# X=Xc[Xc>=np.nanmean(Xc)-np.nanstd(Xc)];Y=Yc[Xc>=np.nanmean(Xc)-np.nanstd(Xc)];
X=Xc[(Xc>=vamp)&(Wc>=vrad)];Y=Yc[(Xc>=vamp)&(Wc>=vrad)];



slope, intercept, r_value, p_value, std_err = stats.linregress(X,Y)
line = slope*(X)+intercept

plt.scatter(X,Y,c=colorx,edgecolors='white',linewidths=1,alpha=0.7)
plt.plot(X,line,color=colorx,alpha=0.5)

# textstr = '\n'.join(('pearson r ={:.2f}'.format(r_value),
#                     'p value ={0}'.format(p_value)))
                    
textstr = ('pearson r ={:.2f}'.format(r_value))


# plt.text(X.min(),Y.max()-0.5,textstr)
plt.title(textstr,fontweight='bold',fontsize=fsz-2)
plt.grid()

plt.ylabel(r'Radius (km)',fontsize=fsz, fontweight='bold')
plt.xlabel('Amplitude (cm)',fontsize=fsz, fontweight='bold')
plt.ylim(40,200)
plt.xlim(0,60)


# ax4 = plt.subplot(gs[0,15:19])

# xxc=Aamp*100
# Xc,Yc=Arad/1e3,Apvel3;
# # X=Xc[Xc>=1];Y=Yc[Xc>=1];
# X=Xc[(xxc>=vamp)&(Wc>=vrad)];Y=Yc[(xxc>=vamp)&(Wc>=vrad)];

# # X=Xc[Xc>=np.nanmean(Xc)-np.nanstd(Xc)];Y=Yc[Xc>=np.nanmean(Xc)-np.nanstd(Xc)];




# slope, intercept, r_value, p_value, std_err = stats.linregress(X,Y)
# line = slope*(X)+intercept

# plt.scatter(X,Y,c=colorx,edgecolors='white',linewidths=1,alpha=0.7)
# plt.plot(X,line,color=colorx,alpha=0.5)

# # textstr = '\n'.join(('pearson r ={:.2f}'.format(r_value),
# #                     'p value ={0}'.format(p_value)))
                    
# textstr = ('pearson r ={:.2f}'.format(r_value))


# # plt.text(X.min(),Y.max()-0.5,textstr)
# plt.title(textstr,fontweight='bold',fontsize=fsz-2)
# plt.grid()

# plt.ylabel(r'Propagating Speed (cm $s^{-1}$)',fontsize=fsz, fontweight='bold')
# plt.xlabel('Radius (km)',fontsize=fsz, fontweight='bold')
# plt.ylim(0,30)
# plt.xlim(vrad,150)

pc.subplots_adjust(top=0.88,
bottom=0.11,
left=0.06,
right=0.97,
hspace=0.2,
wspace=0.2)
plt.tight_layout()

plt.savefig('/home/igor/Dropbox/META_eddydetection_images/Scatter_AntiCyclone.png')
plt.close()


'''
Figura 2 - Histogramas
'''
dbol=True

Ra=Arad/1e3;
Rc=Crad/1e3;

plt.rcParams['font.weight']='bold'
plt.rcParams['font.size']=15
p=plt.figure(figsize=(11.66,  7.81))

gs = gridspec.GridSpec(6,8)
ax1 = plt.subplot(gs[0:3, 0:2])

Hc=-Camp*100
# cH=Hc.copy()[Hc<-1]
cH=Hc.copy()#[(Hc<=-vamp)&(Rc>=40)]

# cH=Hc.copy()[Hc<np.nanmean(Hc)+np.nanstd(Hc)]


Ha=Aamp*100
# AA=Ah[Ha>1]
aH=Ha.copy()#[(Ha>=vamp)&(Ra>=40)]
# aH=Ha.copy()[Ha>np.nanmean(Ha) -np.nanstd(Ha)]


plt.hist(cH,20,edgecolor='k',color='navy',alpha=0.7,density=dbol)
plt.hist(aH,10,edgecolor='k',color='red',alpha=0.7,density=dbol)


if dbol:
    plt.ylabel('Density',fontweight='bold')
else:
    plt.ylabel('Number of eddies',fontweight='bold')

plt.xlabel('Amplitude [cm]',fontweight='bold')


plt.xticks(np.linspace(-30,30,5),np.int_(np.abs(np.linspace(-30,30,5))))

plt.xlim(-30,30)



ax2 = plt.subplot(gs[0:3, 3:5])

Ch=-Cvel*100
# CC=Ch[Hc<-1]
CC=Ch.copy()#[(Hc<=-vamp)&(Rc>=40)]

# CC=Ch[Hc<np.nanmean(Hc)+np.nanstd(Hc)]



Ah=Avel*100
# AA=Ah[Ha>1]
AA=Ah.copy()#[(Ha>=vamp)&(Ra>=40)]

# AA=Ah[Ha>np.nanmean(Ha) -np.nanstd(Ha)]

plt.hist(CC,20,edgecolor='k',color='navy',alpha=0.7,density=dbol)
plt.hist(AA,15,edgecolor='k',color='red',alpha=0.7,density=dbol)
# plt.ylabel('Number of eddies',fontweight='bold')
plt.xlabel(r'Swirling speed [cm s$^{-1}$]',fontweight='bold')

plt.xticks(np.linspace(-60,60,5),np.int_(np.abs(np.linspace(-60,60,5))))
plt.xlim(-60,60)

ax3 = plt.subplot(gs[0:3, 6:8])

Ch=-Ceke
# CC=Ch[Hc<-1]
CC=Ch.copy()#[(Hc<=-vamp)&(Rc>=40)]

# CC=Ch[Hc<np.nanmean(Hc)+np.nanstd(Hc)]


Ah=Aeke
# AA=Ah[Ha>1]
AA=Ah.copy()#[(Ha>=vamp)&(Ra>=40)]

# AA=Ah[Ha>np.nanmean(Ha) -np.nanstd(Ha)]



plt.hist(CC,30,edgecolor='k',color='navy',alpha=0.7,density=dbol)
plt.hist(AA,20,edgecolor='k',color='red',alpha=0.7,density=dbol)
# plt.ylabel('Number of eddies',fontweight='bold')
plt.xlabel('EKE [cm$^{2}$ s$^{-2}$]',fontweight='bold')

plt.xticks(np.linspace(-500,500,3),np.int_(np.abs(np.linspace(-500,500,3))))
plt.xlim(-800,800)


ax4 = plt.subplot(gs[3:6, 0:2])

Ch=-Crad/1e3
# CC=Ch[Hc<-1]
CC=Ch.copy()#[(Hc<=-vamp)&(Rc>=40)]

# CC=Ch[Hc<np.nanmean(Hc)+np.nanstd(Hc)]


Ah=Arad/1e3
# AA=Ah[Ha>1]
AA=Ah.copy()#[(Ha>=vamp)&(Ra>=40)]

# AA=Ah[Ha>np.nanmean(Ha) -np.nanstd(Ha)]

plt.hist(CC,20,edgecolor='k',color='navy',alpha=0.7,density=dbol)
plt.hist(AA,15,edgecolor='k',color='red',alpha=0.7,density=dbol)
# plt.ylabel('Number of eddies',fontweight='bold')
plt.xlabel('Radius [km]',fontweight='bold')

plt.xticks(np.linspace(-200,200,5),np.int_(np.abs(np.linspace(-200,200,5))))
plt.xlim(-200,200)


ax5 = plt.subplot(gs[3:6, 3:5])


# Ch=-Cpvel3
# # CC=Ch[Hc<-1]
# CC=Ch.copy()[(Hc<=-vamp)&(Rc>=40)]

# # CC=Ch[Hc<np.nanmean(Hc)+np.nanstd(Hc)]


# Ah=Apvel3
# # AA=Ah[Ha>1]
# AA=Ah.copy()[(Ha>=vamp)&(Ra>=40)]

# # AA=Ah[Ha>np.nanmean(Ha) -np.nanstd(Ha)]


# plt.hist(CC,10,edgecolor='k',color='navy',alpha=0.7,density=dbol)
# plt.hist(AA,10,edgecolor='k',color='red',alpha=0.7,density=dbol)


# if dbol:
#     plt.ylabel('Density',fontweight='bold')
# else:
#     plt.ylabel('Number of eddies',fontweight='bold')
# plt.xlabel(r'Propagating speed [cm s$^{-1}$]',fontweight='bold')

# plt.xticks(np.linspace(-25,25,5),np.int_(np.abs(np.linspace(-25,25,5))))


# plt.xlim(-25,25)




Ch=-Cdist
# CC=Ch[Hc<-1]
CC=Ch.copy()#[(Hc<=-vamp)&(Rc>=40)]

# CC=Ch[Hc<np.nanmean(Hc)+np.nanstd(Hc)]


Ah=Adist
# AA=Ah[Ha>1]
AA=Ah.copy()#[(Ha>=vamp)&(Ra>=40)]

# AA=Ah[Ha>np.nanmean(Ha) -np.nanstd(Ha)]


plt.hist(CC,15,edgecolor='k',color='navy',alpha=0.7,density=dbol)
plt.hist(AA,15,edgecolor='k',color='red',alpha=0.7,density=dbol)
# plt.ylabel('Number of eddies',fontweight='bold')
plt.xlabel('Distance [km]  $x10^{3}$',fontweight='bold')


plt.xticks(np.linspace(-1500,1500,5),[1.5,1,0,1,1.5],fontsize=12)

plt.xlim(-1500,1500)

ax6 = plt.subplot(gs[3:6, 6:8])



Ch=-Cdur
# CC=Ch[Hc<-1]
CC=Ch.copy()#[(Hc<=-vamp)&(Rc>=40)]

# CC=Ch[Hc<np.nanmean(Hc)+np.nanstd(Hc)]


Ah=Adur
# AA=Ah[Ha>1]
AA=Ah.copy()#[(Ha>=vamp)&(Ra>=40)]

# AA=Ah[Ha>np.nanmean(Ha) -np.nanstd(Ha)]


plt.hist(CC,20,edgecolor='k',color='navy',alpha=0.7,density=dbol)
plt.hist(AA,15,edgecolor='k',color='red',alpha=0.7,density=dbol)
# plt.ylabel('Number of eddies',fontweight='bold')
plt.xlabel('Duration [days]',fontweight='bold')

plt.xticks(np.linspace(-250,250,3),np.int_(np.abs(np.linspace(-250,250,3))))
plt.xlim(-300,300)


plt.tight_layout()
plt.savefig('/home/igor/Dropbox/META_eddydetection_images/Total_histogram.png')
plt.close()



'Figuras 3a,3b,3c - Boxplots'
mostr=[u'Jan', u'Feb', u'Mar', u'Apr',u'May',u'Jun',u'Jul',u'Aug',u'Sep',u'Oct',u'Nov',u'Dec']

'Cyclonic parameters'

Bc=-Camp*100

cmo = pd.DatetimeIndex(np.array(Cdate)[ccond]).month
cfmo= pd.DatetimeIndex(np.array(Cfdate)[ccond]).month
cyea= pd.DatetimeIndex(np.array(Cdate)[ccond]).year
# cyea=cyea[Bc<-vamp];cmo=cmo[Bc<-vamp];cfmo=cfmo[Bc<-vamp];

Cbox=np.zeros([12,len( np.arange(1993, Cdate[-1].year +1 ) )])


for mi,yi in zip(cmo,cyea):
    Cbox[mi-1,yi-1993]=len(cmo[(cmo==mi)&(cyea==yi)])

Cbox=pd.DataFrame(Cbox,columns=np.arange(1993, Cdate[-1].year +1 ),index=mostr)

'Anticyclonic parameters'

Ba=Aamp*100

amo = pd.DatetimeIndex(np.array(Adate)[acond]).month
afmo= pd.DatetimeIndex(np.array(Afdate)[acond]).month
ayea= pd.DatetimeIndex(np.array(Adate)[acond]).year
# amo=amo[Ba>vamp];ayea=ayea[Ba>vamp];afmo=afmo[Ba>vamp]

Abox=np.zeros([12,len( np.arange(1993, Adate[-1].year +1 ) )])


for mi,yi in zip(cmo,cyea):
    Abox[mi-1,yi-1993]=len(amo[(amo==mi)&(ayea==yi)])

Abox=pd.DataFrame(Abox,columns=np.arange(1993, Adate[-1].year +1 ),index=mostr)

mrange=np.arange(0.5,13.5)
yrange=np.arange(cyea[0]-0.5,cyea[-1]+1 +0.5)



plt.rcParams['font.weight']='bold'
plt.rcParams['font.size']=20


'BOX PLOT CYC ANTCYC'
p=plt.figure(figsize=(17.94,  7.46))
gs = gridspec.GridSpec(1,2)
ax1 = plt.subplot(gs[0,0])

cbp=plt.boxplot(Cbox,widths=0.5,patch_artist=True)

for box in cbp['boxes']:
    box.set(color='dimgray',linewidth=2)
    box.set(facecolor ='cornflowerblue')
    # box.set(hatch = '/')

plt.setp(cbp['whiskers'],color='dimgray',linewidth=2)
plt.setp(cbp['fliers'],color='dimgray',linewidth=2)
plt.setp(cbp['medians'],color='dimgray',linewidth=2)
plt.setp(cbp['caps'],color='dimgray',linewidth=2)

# for whisk in cbp
#     whisk.set(color='dimgray',linewidth=2)
#     # fly.set(color='k', linewidth=2,alpha=0.8)


plt.grid()
plt.xticks(ticks=np.arange(1,13),labels=mostr,fontsize=20)
plt.ylim(0,16)


ax2 = plt.subplot(gs[0,1])

abp=plt.boxplot(Abox,widths=0.5,patch_artist=True)

for box in abp['boxes']:
    box.set(color='dimgray',linewidth=2)
    box.set(facecolor = 'lightcoral')
    # box.set(hatch = '/')

plt.setp(abp['whiskers'],color='dimgray',linewidth=2)
plt.setp(abp['fliers'],color='dimgray',linewidth=2)
plt.setp(abp['medians'],color='dimgray',linewidth=2)
plt.setp(abp['caps'],color='dimgray',linewidth=2)


plt.grid()
plt.xticks(ticks=np.arange(1,13),labels=mostr,fontsize=20)
plt.ylim(0,16)

plt.tight_layout()
plt.savefig('/home/igor/Dropbox/META_eddydetection_images/Boxplot_Month.png')
plt.close()



'BOX PLOT SEASONAL'
mosea=['Summer','Autumn','Winter','Spring']

CSazBox=np.zeros([4,28])
ASazBox=np.zeros([4,28])

'Verão'
ASazBox[0,:]=Abox.iloc[:3,:].sum(0) 
CSazBox[0,:]=Cbox.iloc[:3,:].sum(0) 

'Outono'
ASazBox[1,:]=Abox.iloc[3:6,:].sum(0) 
CSazBox[1,:]=Cbox.iloc[3:6,:].sum(0) 

'Inverno'
ASazBox[2,:]=Abox.iloc[6:9,:].sum(0) 
CSazBox[2,:]=Cbox.iloc[6:9,:].sum(0) 

'Primavera'
ASazBox[3,:]=Abox.iloc[9:,:].sum(0) 
CSazBox[3,:]=Cbox.iloc[9:,:].sum(0) 

ASazBox=pd.DataFrame(ASazBox,columns=range(1993,28+1993))
CSazBox=pd.DataFrame(CSazBox,columns=range(1993,28+1993))

ASazBox=ASazBox.iloc[:,:-1]
CSazBox=CSazBox.iloc[:,:-1]



p=plt.figure(figsize=(9.75, 9.47 ))
gs = gridspec.GridSpec(2,1)
ax1 = plt.subplot(gs[0,0])

cbp=plt.boxplot(CSazBox,widths=0.5,patch_artist=True)

for box in cbp['boxes']:
    box.set(color='dimgray',linewidth=2)
    box.set(facecolor ='cornflowerblue')
    # box.set(hatch = '/')

plt.setp(cbp['whiskers'],color='dimgray',linewidth=2)
plt.setp(cbp['fliers'],color='dimgray',linewidth=2)
plt.setp(cbp['medians'],color='dimgray',linewidth=2)
plt.setp(cbp['caps'],color='dimgray',linewidth=2)

plt.text(1,5,'Cyclones')
plt.ylabel('Ocurrences',fontweight='bold')


plt.ylim(2,35)
plt.grid()
plt.xticks(ticks=np.arange(1,5),labels=mosea)
# plt.ylim(0,20)


ax2 = plt.subplot(gs[1,0])

abp=plt.boxplot(ASazBox,widths=0.5,patch_artist=True)

for box in abp['boxes']:
    box.set(color='dimgray',linewidth=2)
    box.set(facecolor = 'lightcoral')
    # box.set(hatch = '/')

plt.setp(abp['whiskers'],color='dimgray',linewidth=2)
plt.setp(abp['fliers'],color='dimgray',linewidth=2)
plt.setp(abp['medians'],color='dimgray',linewidth=2)
plt.setp(abp['caps'],color='dimgray',linewidth=2)

plt.text(1,5,'Anticyclones')
plt.ylabel('Ocurrences',fontweight='bold')

plt.ylim(2,35)

plt.grid()
plt.xticks(ticks=np.arange(1,5),labels=mosea)
# plt.ylim(0,20)

plt.tight_layout()
plt.savefig('/home/igor/Dropbox/META_eddydetection_images/Boxplot_Seasonal.png')
plt.close()



'PDF annual frequency'

p=plt.figure(figsize=(16.51,  6.75 ))
gs = gridspec.GridSpec(1,2)

ax1 = plt.subplot(gs[0, 0])

n=plt.hist(cyea,bins=yrange,color='navy',edgecolor='k',alpha=0.5,density=1)
# plt.plot(np.arange(cyea[0],cyea[-1]+1),n[0],'-',color='navy',alpha=1)
# sns.distplot(cyea,bins=yrange,color='navy',hist_kws={'edgecolor':'k'})
ax1.set_xlabel('Years',fontweight='bold')

ax1.set_xlim(cyea[0]-1,cyea[-1]+1)
# ax1.set_ylim(n[0].min()-5,n[0].max()+5)

ax1.set_xticks(np.arange(cyea[0],cyea[-1]+1)[::5])
ax1.set_xticklabels(np.arange(cyea[0],cyea[-1]+1)[::5])

plt.xlim(1992.5,2019.5)


ax2 = plt.subplot(gs[0, 1])

n=plt.hist(ayea,bins=yrange,color='red',edgecolor='k',alpha=0.5,density=1)
# plt.plot(np.arange(ayea[0],ayea[-1]+1),n[0],'-',color='red',alpha=1)

# sns.distplot(ayea,bins=yrange,color='red',hist_kws={'edgecolor':'k'})
ax2.set_xlabel('Years',fontweight='bold')

ax2.set_xlim(ayea[0]-1,ayea[-1]+1)
# ax2.set_ylim(n[0].min()-5,n[0].max()+5)

ax2.set_xticks(np.arange(ayea[0],ayea[-1]+1)[::5])
ax2.set_xticklabels(np.arange(ayea[0],ayea[-1]+1)[::5])
plt.xlim(1992.5,2019.5)

plt.tight_layout()
plt.savefig('/home/igor/Dropbox/META_eddydetection_images/Total_frequency_anual.png')
plt.close()




'PLotting Formation A+C,C,A Seasonal Fluctuation'
ltmn=-29;ltmx=-14;lnmn=-52;lnmx=-35.

# ltmn=-32;ltmx=-14;lnmn=-52;lnmx=-35.

resol=1.5 #criando mapa de zeros para distribuição
lonf=np.arange(lnmn,lnmx+resol,resol)
latf=np.arange(ltmn,ltmx+resol,resol)

# latf,lonf=np.meshgrid(latf,lonf)
lonf,latf=np.meshgrid(lonf,latf)

'Season setting'
# for ira in range(atrack.max()):
sea=[[1,2,3],
    [4,5,6],
    [7,8,9],
    [10,11,12]]

# i=3
for i in range(4):
    mon=sea[i]


    acmap=np.zeros(lonf.shape);
    ccmap=np.zeros(lonf.shape);


    '''
    FORMATION REGIONS
    '''

    for ira in aunik[acond]:

            # arr=np.mean(anc['effective_radius'][:].data[atrack==ira])/1e3
            # app=np.mean(anc['amplitude'][:].data[atrack==ira])*100
            
            atim=dt0+ timedelta(days=(np.int(np.min(Atime[atrack==ira]))));atim=atim.month;# minimo para pegar o primeiro dia

            # if (arr>=40) & (app>=vamp) & (atim==mon):
            # if (arr>=40) & (app>=vamp) & (atim in mon):
            if (atim in mon):


                alx,aly=alon[atrack==ira][0],alat[atrack==ira][0]
                xind=np.argwhere(alx>=lonf[0,:]).T[0][-1]
                yind=np.argwhere(aly>=latf[:,0]).T[0][-1]
            

                acmap[yind,xind]+=1

    # for irc in range(ctrack.max()):
    for irc in cunik[ccond]:
            # crr=np.mean(cnc['effective_radius'][:].data[ctrack==irc])/1e3
            # cpp=np.mean(cnc['amplitude'][:].data[ctrack==irc])*100
            ctim=dt0+ timedelta(days=(np.int(np.min(Ctime[ctrack==irc]))));ctim=ctim.month; # minimo para pegar o primeiro dia

            if (ctim in mon):
    #        if (crr>=40) & (cpp>=vamp) & (ctim==mon):

                clx,cly=clon[ctrack==irc][0],clat[ctrack==irc][0]
                xind=np.argwhere(clx>=lonf[0,:]).T[0][-1]
                yind=np.argwhere(cly>=latf[:,0]).T[0][-1]

                ccmap[yind,xind]+=1

    p=plt.figure(figsize=(8.06,  6.19))
    lbs=[True, False,False,True]
    C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=-29,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

    C.drawcoastlines(linewidth=1);
    C.fillcontinents(color='#fecf78')
    C.drawstates(linewidth=0.5)
    C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,zorder=100,linewidth=0.1)
    C.drawparallels(np.arange(ltmn,ltmx,4),labels=lbs,zorder=101,linewidth=0.1)

    C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
    # C.contourf(lonf,latf,((acmap+ ccmap)/np.nanmax(acmap + ccmap))*100,cmap='YlGnBu_r',latlon=True,levels=np.linspace(0,60,35))
    # plt.pcolormesh(lonf,latf,((acmap+ ccmap)/np.nanmax(acmap + ccmap))*100,cmap=cm.cm.tempo,vmin=0,vmax=80);
    # plt.pcolormesh(lonf,latf,((acmap+ ccmap)/(np.sum(acmap) + np.sum(ccmap)))*100,cmap=cm.cm.tempo,vmin=0,vmax=5);
    plt.pcolormesh(lonf,latf,((acmap+ ccmap)/(anumber + cnumber))*100,cmap=cm.cm.tempo,vmin=0,vmax=1.5);

    cbar=plt.colorbar()
    cbar.set_label('Formation {0} [%]'.format(mosea[i]),fontweight='bold') 

    C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)

    # plt.title('N = {:d}'.format((np.int(np.sum(ccmap + acmap)))))

    plt.tight_layout()
    # plt.savefig('/home/igor/Documents/AntSimi/figures/Total_formation_season_{0}.png'.format(mosea[i]))
    plt.savefig('/home/igor/Dropbox/META_eddydetection_images/Total_formation_season_{0}_modo1.png'.format(mosea[i]))

    plt.close()



    p=plt.figure(figsize=(8.06,  6.19))
    lbs=[True, False,False,True]
    C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=-29,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

    C.drawcoastlines(linewidth=1);
    C.fillcontinents(color='#fecf78')
    C.drawstates(linewidth=0.5)
    C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,zorder=100,linewidth=0.1)
    C.drawparallels(np.arange(ltmn,ltmx,4),labels=lbs,zorder=101,linewidth=0.1)

    C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
    # C.contourf(lonf,latf,((acmap+ ccmap)/np.nanmax(acmap + ccmap))*100,cmap='YlGnBu_r',latlon=True,levels=np.linspace(0,60,35))
    # plt.pcolormesh(lonf,latf,((acmap+ ccmap)/np.nanmax(acmap + ccmap))*100,cmap=cm.cm.tempo,vmin=0,vmax=80);
    # plt.pcolormesh(lonf,latf,((ccmap)/(np.sum(ccmap)))*100,cmap='Blues',vmin=0,vmax=1);
    plt.pcolormesh(lonf,latf,((ccmap)/(cnumber))*100,cmap='Blues',vmin=0,vmax=2);
    plt.text(lonf.min()+2,latf.max()-2,'{0}'.format(mosea[i]))

    cbar=plt.colorbar()
    cbar.set_label('Formation rate [%]',fontweight='bold') 

    C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)
    # plt.title('N = {:d}'.format((np.int(np.sum(ccmap)))))

    plt.tight_layout()
    # plt.savefig('/home/igor/Documents/AntSimi/figures/Cyclonic_birth_season_{0}.png'.format(mosea[i]))
    plt.savefig('/home/igor/Dropbox/META_eddydetection_images/Cyclonic_birth_season_{0}_modo1.png'.format(mosea[i]))

    plt.close()




    p=plt.figure(figsize=(8.06,  6.19))
    lbs=[True, False,False,True]
    C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=-29,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

    C.drawcoastlines(linewidth=1);
    C.fillcontinents(color='#fecf78')
    C.drawstates(linewidth=0.5)
    C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,zorder=100,linewidth=0.1)
    C.drawparallels(np.arange(ltmn,ltmx,4),labels=lbs,zorder=101,linewidth=0.1)

    C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
    # C.contourf(lonf,latf,((acmap+ ccmap)/np.nanmax(acmap + ccmap))*100,cmap='YlGnBu_r',latlon=True,levels=np.linspace(0,60,35))
    # plt.pcolormesh(lonf,latf,((acmap+ ccmap)/np.nanmax(acmap + ccmap))*100,cmap=cm.cm.tempo,vmin=0,vmax=80);
    # plt.pcolormesh(lonf,latf,((acmap)/(np.sum(acmap)))*100,cmap='Reds',vmin=0,vmax=1);
    plt.pcolormesh(lonf,latf,((acmap)/(anumber))*100,cmap='Reds',vmin=0,vmax=2);
    plt.text(lonf.min()+2,latf.max()-2,'{0}'.format(mosea[i]))
    cbar=plt.colorbar()
    cbar.set_label('Formation rate [%]',fontweight='bold') 

    C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)
    # plt.title('N = {:d}'.format((np.int(np.sum(acmap)))))


    plt.tight_layout()
    # plt.savefig('/home/igor/Documents/AntSimi/figures/Anticyclonic_birth_season_{0}.png'.format(mosea[i]))
    plt.savefig('/home/igor/Dropbox/META_eddydetection_images/Anticyclonic_birth_season_{0}_modo1.png'.format(mosea[i]))
    plt.close()







'Figures 4 - Heatmaps of Dissipation, Formation, etc'

ltmn=-32;ltmx=-14;lnmn=-52;lnmx=-35.

resol=1.5
lonf=np.arange(lnmn,lnmx+resol,resol)
latf=np.arange(ltmn,ltmx+resol,resol)

# latf,lonf=np.meshgrid(latf,lonf)
lonf,latf=np.meshgrid(lonf,latf)

acmap=np.zeros(lonf.shape);
ccmap=np.zeros(lonf.shape);



'''
FORMATION REGIONS
'''
# for ira in range(atrack.max()):
for ira in aunik[acond]:

        arr=np.mean(anc['effective_radius'][:].data[atrack==ira])/1e3
        app=np.mean(anc['amplitude'][:].data[atrack==ira])*100

        if (arr>=40) & (app>=vamp):

            alx,aly=alon[atrack==ira][0],alat[atrack==ira][0]
            xind=np.argwhere(alx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(aly>=latf[:,0]).T[0][-1]


            acmap[yind,xind]+=1

# for irc in range(ctrack.max()):
for irc in cunik[ccond]:
        crr=np.mean(cnc['effective_radius'][:].data[ctrack==irc])/1e3
        cpp=np.mean(cnc['amplitude'][:].data[ctrack==irc])*100

        if (crr>=40) & (cpp>=vamp):

            clx,cly=clon[ctrack==irc][0],clat[ctrack==irc][0]
            xind=np.argwhere(clx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(cly>=latf[:,0]).T[0][-1]

            ccmap[yind,xind]+=1



p=plt.figure(figsize=(8.06,  6.19))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=-29,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,4),labels=lbs,zorder=101,linewidth=0.1)

C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
# C.contourf(lonf,latf,((acmap+ ccmap)/np.nanmax(acmap + ccmap))*100,cmap='YlGnBu_r',latlon=True,levels=np.linspace(0,60,35))
# plt.pcolormesh(lonf,latf,((acmap+ ccmap)/np.nanmax(acmap + ccmap))*100,cmap=cm.cm.tempo,vmin=0,vmax=80);
plt.pcolormesh(lonf,latf,((acmap+ ccmap)/(anumber + cnumber))*100,cmap=cm.cm.tempo,vmin=0,vmax=6);

cbar=plt.colorbar()
cbar.set_label('Formation [%]',fontweight='bold') 

C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)


plt.tight_layout()
# plt.savefig('/home/igor/Documents/AntSimi/figures/Total_formation.png')
plt.savefig('/home/igor/Dropbox/META_eddydetection_images/Total_formation_modo2.png')

plt.close()



'''
DISSIPATION REGIONS
'''
acdmap=np.zeros(lonf.shape);
ccdmap=np.zeros(lonf.shape);


# for ira in range(atrack.max()):
for ira in aunik[acond]:

        arr=np.mean(anc['effective_radius'][:].data[atrack==ira])/1e3
        app=np.mean(anc['amplitude'][:].data[atrack==ira])*100

        if (arr>=40) & (app>=vamp):

            alx,aly=alon[atrack==ira][-1],alat[atrack==ira][-1]
            xind=np.argwhere(alx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(aly>=latf[:,0]).T[0][-1]

            acdmap[yind,xind]+=1

# for irc in range(ctrack.max()):
for irc in cunik[ccond]:
        crr=np.mean(cnc['effective_radius'][:].data[ctrack==irc])/1e3
        cpp=np.mean(cnc['amplitude'][:].data[ctrack==irc])*100

        if (crr>=40) & (cpp>=vamp):
            clx,cly=clon[ctrack==irc][-1],clat[ctrack==irc][-1]
            xind=np.argwhere(clx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(cly>=latf[:,0]).T[0][-1]

            ccdmap[yind,xind]+=1


p=plt.figure(figsize=(8.06,  6.19))
lbs=[True, False,False,True]

C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=-29,urcrnrlat=ltmx-1,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,4),labels=lbs,zorder=101,linewidth=0.1)

C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
# C.contourf(lonf,latf,(acdmap + ccdmap)/np.nanmax(acdmap + ccdmap)*100,cmap='YlGnBu_r',latlon=True,levels=np.linspace(0,60,35))
# plt.pcolormesh(lonf,latf,((acdmap+ ccdmap)/np.nanmax(acdmap + ccdmap))*100,cmap=cm.cm.tempo,vmin=0,vmax=80);
plt.pcolormesh(lonf,latf,((acdmap+ ccdmap)/(anumber + cnumber))*100,cmap=cm.cm.tempo,vmin=0,vmax=5);

cbar=plt.colorbar()
cbar.set_label('Dissipation [%]',fontweight='bold') 
C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)
plt.tight_layout()
plt.savefig('/home/igor/Dropbox/META_eddydetection_images/Total_dissipation_modo2.png')
plt.close()


# '''
# POLARITY BIRTH MAP
# '''

# 'BOTH'
# p=plt.figure(figsize=(8.06,  6.19))
# lbs=[True, False,False,True]
# C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=-29,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

# C.drawcoastlines(linewidth=1);
# C.fillcontinents(color='#fecf78')
# C.drawstates(linewidth=0.5)
# C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,zorder=100,linewidth=0.1)
# C.drawparallels(np.arange(ltmn,ltmx,4),labels=lbs,zorder=101,linewidth=0.1)

# C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,200],linewidths=0.5,latlon=True)
# # C.contourf(lonf,latf,(acmap - ccmap)/np.nanmax(abs(acmap - ccmap))*100,cmap='seismic',latlon=True,levels=np.linspace(-100,100,35))
# C.contourf(lonf,latf,((acmap - ccmap)/np.sum(acmap-ccmap))*100,cmap='seismic',latlon=True,levels=np.linspace(-100,100,35))

# # plt.pcolormesh(lonf,latf,((acmap - ccmap)/np.nanmax(acmap - ccmap))*100,cmap='seismic',vmin=-100,vmax=100);

# cbar=plt.colorbar(ticks=np.arange(-100,120,20))
# cbar.set_label('Polarity Formation [%]',fontweight='bold')

# C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)
# plt.tight_layout()
# plt.savefig('/home/igor/Documents/AntSimi/figures/Total_polarity_birth_modo2.png')
# plt.close()


'CYCLONES'
p=plt.figure(figsize=(8.06,  6.19))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=-29,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,4),labels=lbs,zorder=101,linewidth=0.1)

C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
# C.contourf(lonf,latf,(acmap - ccmap)/np.nanmax(abs(acmap - ccmap))*100,cmap='seismic',latlon=True,levels=np.linspace(-100,100,35))
# plt.pcolormesh(lonf,latf,((ccmap)/np.nanmax(ccmap))*100,cmap='Blues',vmin=0,vmax=80);
plt.pcolormesh(lonf,latf,((ccmap)/cnumber)*100,cmap='Blues',vmin=0,vmax=6);


cbar=plt.colorbar()
cbar.set_label(' Cyclone Formation [%]',fontweight='bold')

C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)
plt.tight_layout()
plt.savefig('/home/igor/Dropbox/META_eddydetection_images/cyclonic_birth.png')
plt.close()


'ANTICYCLONES'
p=plt.figure(figsize=(8.06,  6.19))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=-29,urcrnrlat=ltmx,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,4),labels=lbs,zorder=101,linewidth=0.1)

C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
# C.contourf(lonf,latf,(acmap - ccmap)/np.nanmax(abs(acmap - ccmap))*100,cmap='seismic',latlon=True,levels=np.linspace(-100,100,35))
plt.pcolormesh(lonf,latf,((acmap)/anumber)*100,cmap='Reds',vmin=0,vmax=6);

cbar=plt.colorbar()
cbar.set_label(' Anticyclone Formation [%]',fontweight='bold')

C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)
plt.tight_layout()
plt.savefig('/home/igor/Dropbox/META_eddydetection_images/anticyclonic_birth.png')
plt.close()



'''
SWRILING SPEED MAP
'''

velmap=np.zeros(lonf.shape);
adummy=np.zeros(lonf.shape);
cdummy=np.zeros(lonf.shape);


# for ira in range(atrack.max()):
for ira in aunik[acond]:
        arr=np.mean(anc['effective_radius'][:].data[atrack==ira])/1e3
        app=np.mean(anc['amplitude'][:].data[atrack==ira])*100

        if (arr>=40) & (app>=vamp):

            alx,aly=np.mean(alon[atrack==ira]),np.mean(alat[atrack==ira])

            xind=np.argwhere(alx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(aly>=latf[:,0]).T[0][-1]

            adummy[yind,xind]+=1
            velmap[yind,xind]+=Avel[ira]*100


# for irc in range(ctrack.max()):
for irc in cunik[ccond]:
        crr=np.mean(cnc['effective_radius'][:].data[ctrack==irc])/1e3
        cpp=np.mean(cnc['amplitude'][:].data[ctrack==irc])*100

        if (crr>=40) & (cpp>=vamp):

            clx,cly=np.mean(clon[ctrack==irc]),np.mean(clat[ctrack==irc])
            xind=np.argwhere(clx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(cly>=latf[:,0]).T[0][-1]

            cdummy[yind,xind]+=1
            velmap[yind,xind]+=Cvel[irc]*100

dummy=adummy+cdummy

Velp=np.nanmean(((velmap)/dummy),axis=1)
Velpstd=np.nanstd(((velmap)/dummy),axis=1)

ps=plt.figure(figsize=(5.18, 7.43))
plt.plot(Velp,latf[:,0],color='green')
plt.fill_betweenx(x1=Velp-Velpstd,x2=Velp+Velpstd,y=latf[:,0],color='teal',alpha=0.5) 
plt.grid()
plt.xlabel('Swirling Speed [cm s$^{-1}$]',fontweight='bold')
plt.ylabel('Latitude',fontweight='bold')
plt.ylim(-32,-14)

# 'Latitudes'
# plt.plot([Velpstd.min(),Velpstd.max()],[-20,-20],'k--')



plt.tight_layout()
plt.savefig('/home/igor/Dropbox/META_eddydetection_images/Profile_Total_velocity.png')
plt.close()





scafac=15
'Scatter MAP'
p=plt.figure(figsize=(10.65,  7.43))
C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=-29,urcrnrlat=ltmx-1,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,4),labels=lbs,linewidth=0.1)

# C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)

# plt.pcolormesh(lonf,latf,((velmap)/(dummy)),cmap=cm.cm.tempo,vmin=9,vmax=26);
axx=plt.scatter(lonf,latf,c=((velmap)/(dummy)),s=dummy*scafac,cmap=cm.cm.tempo,edgecolors='k',vmax=25)


handles, labels = axx.legend_elements(prop="sizes", alpha=0.6)
plt.legend(handles[::3],['15','65','115'], loc="lower right", title="Samples") 

cbar=plt.colorbar()
cbar.set_label('Swirling Speed [cm s$^{-1}$]',fontweight='bold') 

plt.tight_layout()



plt.savefig('/home/igor/Dropbox/META_eddydetection_images/SCattermap_Total_velocity.png')
plt.close()



dummy[dummy<5]=0
velmap[dummy<5]=0

p=plt.figure(figsize=(10.65,  7.43))
C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=-29,urcrnrlat=ltmx-1,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,4),labels=lbs,zorder=101,linewidth=0.1)

C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
# C.contourf(lonf,latf,(radmap)/(adummy+cdummy),cmap='rainbow',extend='both',latlon=True,levels=np.linspace(40,130,35))
# plt.pcolormesh(lonf,latf,(radmap)/(adummy+cdummy),cmap='rainbow',vmin=40,vmax=130);
plt.pcolormesh(lonf,latf,((velmap)/(dummy)),cmap=cm.cm.tempo,vmin=9,vmax=26);

cbar=plt.colorbar()
cbar.set_label('Swirling Speed [cm s$^{-1}$]',fontweight='bold') 
C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)
plt.tight_layout()

plt.savefig('/home/igor/Dropbox/META_eddydetection_images/Total_velocity.png')
plt.close()



'''
RADIUS MAP
'''

radmap=np.zeros(lonf.shape);
adummy=np.zeros(lonf.shape);
cdummy=np.zeros(lonf.shape);


# for ira in range(atrack.max()):
for ira in aunik[acond]:
        arr=np.mean(anc['effective_radius'][:].data[atrack==ira])/1e3
        app=np.mean(anc['amplitude'][:].data[atrack==ira])*100

        if (arr>=40) & (app>=vamp):

            alx,aly=np.mean(alon[atrack==ira]),np.mean(alat[atrack==ira])

            xind=np.argwhere(alx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(aly>=latf[:,0]).T[0][-1]

            adummy[yind,xind]+=1
            radmap[yind,xind]+=Arad[ira]/1e3


# for irc in range(ctrack.max()):
for irc in cunik[ccond]:
        crr=np.mean(cnc['effective_radius'][:].data[ctrack==irc])/1e3
        cpp=np.mean(cnc['amplitude'][:].data[ctrack==irc])*100

        if (crr>=40) & (cpp>=vamp):

            clx,cly=np.mean(clon[ctrack==irc]),np.mean(clat[ctrack==irc])
            xind=np.argwhere(clx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(cly>=latf[:,0]).T[0][-1]

            cdummy[yind,xind]+=1
            radmap[yind,xind]+=Crad[irc]/1e3


dummy=adummy+cdummy



Radp=np.nanmean(((radmap)/dummy),axis=1)
Radpstd=np.nanstd(((radmap)/dummy),axis=1)

ps=plt.figure(figsize=(5.18, 7.43))
plt.plot(Radp,latf[:,0],color='green')
plt.fill_betweenx(x1=Radp-Radpstd,x2=Radp+Radpstd,y=latf[:,0],color='teal',alpha=0.5) 
plt.grid()
plt.xlabel('Radii [km]',fontweight='bold')
plt.ylabel('Latitude',fontweight='bold')
plt.ylim(-32,-14)

# 'Latitudes'
# plt.plot([Velpstd.min(),Velpstd.max()],[-20,-20],'k--')



plt.tight_layout()
plt.savefig('/home/igor/Dropbox/META_eddydetection_images/Profile_Total_radii.png')
plt.close()




'Scatter MAP'
p=plt.figure(figsize=(10.65,  7.43))
C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=-29,urcrnrlat=ltmx-1,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,4),labels=lbs,linewidth=0.1)

# C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)

# plt.pcolormesh(lonf,latf,((velmap)/(dummy)),cmap=cm.cm.tempo,vmin=9,vmax=26);
plt.scatter(lonf,latf,c=((radmap)/(dummy)),s=dummy*scafac,vmin=50,vmax=75,cmap=cm.cm.tempo,edgecolors='k')

cbar=plt.colorbar()
cbar.set_label('Radius [km]',fontweight='bold') 

plt.tight_layout()
plt.savefig('/home/igor/Dropbox/META_eddydetection_images/SCattermap_Total_radius.png')
plt.close()




dummy[dummy<5]=0
radmap[dummy<5]=0

p=plt.figure(figsize=(10.65,  7.43))
C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=-29,urcrnrlat=ltmx-1,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,4),labels=lbs,zorder=101,linewidth=0.1)

C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
# C.contourf(lonf,latf,(radmap)/(adummy+cdummy),cmap='rainbow',extend='both',latlon=True,levels=np.linspace(40,130,35))
# plt.pcolormesh(lonf,latf,(radmap)/(adummy+cdummy),cmap='rainbow',vmin=40,vmax=130);
plt.pcolormesh(lonf,latf,((radmap)/(dummy)),cmap=cm.cm.tempo,vmin=40,vmax=90);

cbar=plt.colorbar()
cbar.set_label('Radius [km]',fontweight='bold') 
C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)
plt.tight_layout()

plt.savefig('/home/igor/Dropbox/META_eddydetection_images/Total_radius.png')
plt.close()



'''
AMP MAP
'''

ampmap=np.zeros(lonf.shape);
adummy=np.zeros(lonf.shape);
cdummy=np.zeros(lonf.shape);


# for ira in range(atrack.max()):
for ira in aunik[acond]:
        arr=np.mean(anc['effective_radius'][:].data[atrack==ira])/1e3
        app=np.mean(anc['amplitude'][:].data[atrack==ira])*100

        if (arr>=40) & (app>=vamp):

            alx,aly=np.mean(alon[atrack==ira]),np.mean(alat[atrack==ira])

            xind=np.argwhere(alx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(aly>=latf[:,0]).T[0][-1]

            adummy[yind,xind]+=1
            ampmap[yind,xind]+=Aamp[ira]*100


# for irc in range(ctrack.max()):
for irc in cunik[ccond]:
        crr=np.mean(cnc['effective_radius'][:].data[ctrack==irc])/1e3
        cpp=np.mean(cnc['amplitude'][:].data[ctrack==irc])*100

        if (crr>=40) & (cpp>=vamp):

            clx,cly=np.mean(clon[ctrack==irc]),np.mean(clat[ctrack==irc])
            xind=np.argwhere(clx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(cly>=latf[:,0]).T[0][-1]

            cdummy[yind,xind]+=1
            ampmap[yind,xind]+=Camp[irc]*100




dummy=adummy+cdummy



Ampp=np.nanmean(((ampmap)/dummy),axis=1)
Amppstd=np.nanstd(((ampmap)/dummy),axis=1)

ps=plt.figure(figsize=(5.18, 7.43))
plt.plot(Ampp,latf[:,0],color='green')
plt.fill_betweenx(x1=Ampp-Amppstd,x2=Ampp+Amppstd,y=latf[:,0],color='teal',alpha=0.5) 
plt.grid()
plt.xlabel('Amplitude [cm]',fontweight='bold')
plt.ylabel('Latitude',fontweight='bold')
plt.ylim(-32,-14)

# 'Latitudes'
# plt.plot([Velpstd.min(),Velpstd.max()],[-20,-20],'k--')



plt.tight_layout()
plt.savefig('/home/igor/Dropbox/META_eddydetection_images/Profile_Total_amplitude.png')
plt.close()



'Scatter MAP'
p=plt.figure(figsize=(10.65,  7.43))
C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=-29,urcrnrlat=ltmx-1,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,4),labels=lbs,linewidth=0.1)

# C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)

# plt.pcolormesh(lonf,latf,((velmap)/(dummy)),cmap=cm.cm.tempo,vmin=9,vmax=26);
plt.scatter(lonf,latf,c=((ampmap)/(dummy)),s=dummy*scafac,cmap=cm.cm.tempo,edgecolors='k',vmin=1,vmax=6)

cbar=plt.colorbar()
cbar.set_label('Amplitude [cm]',fontweight='bold') 

plt.tight_layout()
plt.savefig('/home/igor/Dropbox/META_eddydetection_images/SCattermap_Total_amplitude.png')
plt.close()


dummy[dummy<5]=0
ampmap[dummy<5]=0



p=plt.figure(figsize=(10.65,  7.43))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=-29,urcrnrlat=ltmx-1,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,4),labels=lbs,zorder=101,linewidth=0.1)

C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
# C.contourf(lonf,latf,(ampmap)/(adummy+cdummy),cmap='rainbow',extend='both',latlon=True,levels=np.linspace(vamp,10,35))
plt.pcolormesh(lonf,latf,(ampmap)/(dummy),cmap=cm.cm.tempo,vmin=1,vmax=7);

cbar=plt.colorbar()
cbar.set_label('Amplitude [cm]',fontweight='bold') 
C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)
plt.tight_layout()

plt.savefig('/home/igor/Dropbox/META_eddydetection_images/Total_amplitude.png')
plt.close()



'''
NON-LINEARITY
'''


nlmap=np.zeros(lonf.shape);
adummy=np.zeros(lonf.shape);
cdummy=np.zeros(lonf.shape);


# for ira in range(atrack.max()):
i=0
for ira in aunik[acond]:
        arr=np.mean(anc['effective_radius'][:].data[atrack==ira])/1e3
        app=np.mean(anc['amplitude'][:].data[atrack==ira])*100

        if (arr>=40) & (app>=vamp):

            alx,aly=np.mean(alon[atrack==ira]),np.mean(alat[atrack==ira])

            xind=np.argwhere(alx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(aly>=latf[:,0]).T[0][-1]

            adummy[yind,xind]+=1
            # nlmap[yind,xind]+= (np.nanmax(anc['speed_average'][:][atrack==ira])  /  (Apvel[i]))
            nlmap[yind,xind]+= (np.nanmax(anc['speed_average'][:][atrack==ira])  /  (Apvel3[i]/100))

            i+=1

i=0
# for irc in range(ctrack.max()):
for irc in cunik[ccond]:
    crr=np.mean(cnc['effective_radius'][:].data[ctrack==irc])/1e3
    cpp=np.mean(cnc['amplitude'][:].data[ctrack==irc])*100

    if (crr>=40) & (cpp>=vamp):
        clx,cly=np.mean(clon[ctrack==irc]),np.mean(clat[ctrack==irc])

        xind=np.argwhere(clx>=lonf[0,:]).T[0][-1]
        yind=np.argwhere(cly>=latf[:,0]).T[0][-1]

        cdummy[yind,xind]+=1
        # nlmap[yind,xind]+=(np.nanmax(cnc['speed_average'][:][ctrack==irc])/(Cpvel[i]))
        nlmap[yind,xind]+=(np.nanmax(cnc['speed_average'][:][ctrack==irc])/(Cpvel3[i]/100))

        
        i+=1



dummy=adummy+cdummy



Nlp=np.nanmean(((nlmap)/dummy),axis=1)
Nlpstd=np.nanstd(((nlmap)/dummy),axis=1)

ps=plt.figure(figsize=(5.18, 7.43))
plt.plot(Nlp,latf[:,0],color='green')
plt.fill_betweenx(x1=Nlp-Nlpstd,x2=Nlp+Nlpstd,y=latf[:,0],color='teal',alpha=0.5) 
plt.grid()
plt.xlabel('Non-linear Advective Parameter',fontweight='bold')
plt.ylabel('Latitude',fontweight='bold')
plt.ylim(-32,-14)

# 'Latitudes'
# plt.plot([Velpstd.min(),Velpstd.max()],[-20,-20],'k--')



plt.tight_layout()
plt.savefig('/home/igor/Dropbox/META_eddydetection_images/Profile_Total_nonlinear.png')
plt.close()




'Scatter MAP'
p=plt.figure(figsize=(10.65,  7.43))
C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=-29,urcrnrlat=ltmx-1,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,4),labels=lbs,linewidth=0.1)

# C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)

# plt.pcolormesh(lonf,latf,((velmap)/(dummy)),cmap=cm.cm.tempo,vmin=9,vmax=26);
plt.scatter(lonf,latf,c=((nlmap)/(dummy)),s=dummy*scafac,cmap=cm.cm.tempo,edgecolors='k',vmin=2,vmax=15)

cbar=plt.colorbar()
cbar.set_label('Advective Nonlinearity Parameter',fontweight='bold') 

plt.tight_layout()
plt.savefig('/home/igor/Dropbox/META_eddydetection_images/SCattermap_Total_nonlinear.png')
plt.close()


dummy[dummy<5]=0
nlmap[dummy<5]=0


p=plt.figure(figsize=(10.65,  7.43))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=-29,urcrnrlat=ltmx-1,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,4),labels=lbs,zorder=101,linewidth=0.1)

C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
# C.contourf(lonf,latf,(nlmap)/(adummy+cdummy),cmap='rainbow',extend='both',latlon=True,levels=np.linspace(0,6,35))
plt.pcolormesh(lonf,latf,(nlmap)/(dummy),cmap=cm.cm.tempo,vmin=0.5,vmax=5.4);

cbar=plt.colorbar()
cbar.set_label('Advective Nonlinearity Parameter',fontweight='bold') 
C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)
plt.tight_layout()
plt.savefig('/home/igor/Dropbox/META_eddydetection_images/Total_nonlinear.png')
plt.close()






'''
DURATION
'''


durmap=np.zeros(lonf.shape);
adummy=np.zeros(lonf.shape);
cdummy=np.zeros(lonf.shape);


# for ira in range(atrack.max()):
i=0
for ira in aunik[acond]:
        arr=np.mean(anc['effective_radius'][:].data[atrack==ira])/1e3
        app=np.mean(anc['amplitude'][:].data[atrack==ira])*100

        if (arr>=40) & (app>=vamp):

            alx,aly=np.mean(alon[atrack==ira]),np.mean(alat[atrack==ira])

            xind=np.argwhere(alx>=lonf[0,:]).T[0][-1]
            yind=np.argwhere(aly>=latf[:,0]).T[0][-1]

            adummy[yind,xind]+=1
            # nlmap[yind,xind]+= (np.nanmax(anc['speed_average'][:][atrack==ira])  /  (Apvel[i]))
            durmap[yind,xind]+= len(atrack[atrack==ira])

            i+=1

i=0
# for irc in range(ctrack.max()):
for irc in cunik[ccond]:
    crr=np.mean(cnc['effective_radius'][:].data[ctrack==irc])/1e3
    cpp=np.mean(cnc['amplitude'][:].data[ctrack==irc])*100

    if (crr>=40) & (cpp>=vamp):
        clx,cly=np.mean(clon[ctrack==irc]),np.mean(clat[ctrack==irc])

        xind=np.argwhere(clx>=lonf[0,:]).T[0][-1]
        yind=np.argwhere(cly>=latf[:,0]).T[0][-1]

        cdummy[yind,xind]+=1
        # nlmap[yind,xind]+=(np.nanmax(cnc['speed_average'][:][ctrack==irc])/(Cpvel[i]))
        durmap[yind,xind]+=len(ctrack[ctrack==irc])
        
        i+=1



dummy=adummy+cdummy





'Scatter MAP'
p=plt.figure(figsize=(10.65,  7.43))
C = Basemap(llcrnrlon=lnmn+1,llcrnrlat=-29,urcrnrlat=ltmx-1,urcrnrlon=lnmx,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn,lnmx,4),labels=lbs,linewidth=0.1)
C.drawparallels(np.arange(ltmn,ltmx,4),labels=lbs,linewidth=0.1)

# C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)

# plt.pcolormesh(lonf,latf,((velmap)/(dummy)),cmap=cm.cm.tempo,vmin=9,vmax=26);
plt.scatter(lonf,latf,c=((durmap)/(dummy)),s=dummy*scafac,cmap=cm.cm.tempo,edgecolors='k',vmin=15,vmax=90)

cbar=plt.colorbar()
cbar.set_label('Duration [day]',fontweight='bold') 

plt.tight_layout()
plt.savefig('/home/igor/Dropbox/META_eddydetection_images/SCattermap_Total_duration.png')
plt.close()

