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
import matplotlib.patches as patches


'Climatology'

var=Dataset('/home/igor/Documents/AntSimi/altdata/clim/adt_climatology.nc')
dt0=datetime(1950,1,1)
dt=pd.DatetimeIndex([dt0 + timedelta(days=np.int(dia)) for dia in var['time']])
'''
lonc,latc=var['longitude'][:]-360,var['latitude'][:];
slac=var['sla'][:]
adtc=var['adt'][:]
'''

'Topography'
topo=Dataset('/home/igor/Documents/Extras/Sumida/HYCOM/bat/etopo1_bedrock.nc')
xb=topo['lon'][:];yb=topo['lat'][:];zb=topo['Band1'][:];
xxb,yyb=np.meshgrid(xb,yb)


'Track data'
anc=Dataset('/home/igor/Documents/AntSimi_adt/tracking/Anticyclonic.nc')
cnc=Dataset('/home/igor/Documents/AntSimi_adt/tracking/Cyclonic.nc')



edatabase=[[-16,-13,-39.5,-37,anc,1,'lightsalmon','IE'],
            [-18,-16,-39,-36.5,anc,1,'chocolate','RCE'],
            [-20,-18,-38.5,-36,anc,1,'crimson','AE'],
            [-21.4,-19.25,-40.725,-37.750,cnc,-1,'navy','VE'],
            [-22.5,-20.5,-40.725,-38.250,anc,1,'darkorange','SVE'],
            [-24,-21.5,-41,-39,cnc,-1,'cornflowerblue','CSTE'],
            [-25,-23,-43,-41,cnc,-1,'mediumseagreen','CFE'],]

Table=pd.DataFrame([])


'MAP of REGIONS'
'CYCLONES'

plt.rcParams['font.size']=13
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.weight']='bold'

ltmns=-28;ltmxs=-11;lnmns=-46;lnmxs=-35.

fig =plt.figure(figsize=(6.74, 9.22))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmns,llcrnrlat=ltmns,urcrnrlat=ltmxs,urcrnrlon=lnmxs,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmns,lnmxs,2),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmns,ltmxs,2),labels=lbs,zorder=101,linewidth=0.1)



for ity in range(len(edatabase)):
    ltmnw,ltmxw,lnmnw,lnmxw,_,pol,colorw,strw=edatabase[ity]
    if pol==-1:
        plt.fill_between([lnmnw,lnmxw],ltmnw,ltmxw,color=colorw,label=strw,alpha=0.7) 

C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)

plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi_adt/figures/Region_MAP_cyc.png')
plt.close()


'ANTICYCLONES'

plt.rcParams['font.size']=13
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.weight']='bold'

ltmns=-28;ltmxs=-11;lnmns=-46;lnmxs=-35.

fig =plt.figure(figsize=(6.74, 9.22))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmns,llcrnrlat=ltmns,urcrnrlat=ltmxs,urcrnrlon=lnmxs,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmns,lnmxs,2),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmns,ltmxs,2),labels=lbs,zorder=101,linewidth=0.1)



for ity in range(len(edatabase)):
    ltmnw,ltmxw,lnmnw,lnmxw,_,pol,colorw,strw=edatabase[ity]
    if pol==1:
        plt.fill_between([lnmnw,lnmxw],ltmnw,ltmxw,color=colorw,label=strw,alpha=0.7) 

C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=0.5,latlon=True)
C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,np.inf],latlon=True)

plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi_adt/figures/Region_MAP_anti.png')
plt.close()





for itt in range(len(edatabase)):
    ltmn,ltmx,lnmn,lnmx,nnc,polx,colorx,streddy=edatabase[itt]

    xtrack=nnc['track'][:]
    xlon,xlat=nnc['longitude'][:]-360,nnc['latitude'][:]


    ntrack=np.empty(xtrack.shape)*np.nan; #newtrack

    for ix in range(xtrack.max()):
        xts=xlon[xtrack==ix][0];
        yts=xlat[xtrack==ix][0];

        cond=(xts>=lnmn)&(xts<=lnmx)&(yts<=ltmx)&(yts>=ltmn)    

        if cond:
            ntrack[xtrack==ix]=xtrack[xtrack==ix]


    uni=np.unique(ntrack);uni=np.int_(uni[~np.isnan(uni)]);


    Xtime=nnc['time'][:][~np.isnan(ntrack)]
    xtime= [datetime(1950,1,1)+ timedelta(days=np.int(tem)) for tem in Xtime] 
    Xamp=pd.Series(nnc['amplitude'][:]);Xamp=Xamp.groupby(ntrack).mean();
    Xvel=pd.Series(nnc['speed_average'][:]);Xvel=Xvel.groupby(ntrack).mean();
    Xrad=pd.Series(nnc['effective_radius'][:]);Xrad=Xrad.groupby(ntrack).mean();
    Xeke=0.5*((Xvel*100)**2)

    Xdist=[];Xdur=[];Xpvel=[];Xpvel2=[];
    Xdate=[];Xfdate=[];


    for ui in uni:

        dur=len(ntrack[ntrack==ui])

        distx=sw.dist(lon=xlon[ntrack==ui],lat=xlat[ntrack==ui])[0]
        disty=sw.dist(lon=[xlon[ntrack==ui][0],xlon[ntrack==ui][-1]],
                lat=[xlat[ntrack==ui][0],xlat[ntrack==ui][-1]])[0]

        k=np.insert(distx,0,0)
        dist=np.cumsum(k)

        pvel=(dist[-1]*1e3)/(dur*86400)
        # pvel2=np.nanmean((distx*1e3)/(86400))
        pvel2=(disty[0]*1e3/(dur*86400))


        time=nnc['time'][ntrack==ui];

        date=dt0+ timedelta(days=np.int(time.mean()));
        fdate=dt0 + timedelta(days=np.int(time[0]));

        Xdist.append(disty[0]);
        # Xdist.append(dist[-1]);
        Xdur.append(dur);
        Xpvel.append(pvel);
        Xpvel2.append(pvel2);
        Xdate.append(date);
        Xfdate.append(fdate);

    Xdist=np.array(Xdist);
    Xdur=np.array(Xdur);
    Xpvel=np.array(Xpvel);
    Xpvel2=np.array(Xpvel2);

    vamp=1
    vrad=40
    Wc=Xrad/1e3


    xcond=(Xamp>=vamp/100)&(Xrad>=vrad*1e3)

    'Table'
    ttt=pd.Series({'amp':np.nanmean(Xamp[xcond]),'amp_min':np.nanmin(Xamp[xcond]),'amp_max':np.nanmax(Xamp[xcond]),'amp_std':np.nanstd(Xamp[xcond]),
                    'vel':np.nanmean(Xvel[xcond]),'vel_min':np.nanmin(Xvel[xcond]),'vel_max':np.nanmax(Xvel[xcond]),'vel_std':np.nanstd(Xvel[xcond]),
                    'rad':np.nanmean(Xrad[xcond]),'rad_min':np.nanmin(Xrad[xcond]),'rad_max':np.nanmax(Xrad[xcond]),'rad_std':np.nanstd(Xrad[xcond]),
                    'eke':np.nanmean(Xeke[xcond]),'eke_min':np.nanmin(Xeke[xcond]),'eke_max':np.nanmax(Xeke[xcond]),'eke_std':np.nanstd(Xeke[xcond]),
                    'dist':np.nanmean(Xdist[xcond]),'dist_min':np.nanmin(Xdist[xcond]),'dist_max':np.nanmax(Xdist[xcond]),'dist_std':np.nanstd(Xdist[xcond]),
                    'dur':np.nanmean(Xdur[xcond]),'dur_min':np.nanmin(Xdur[xcond]),'dur_max':np.nanmax(Xdur[xcond]),'dur_std':np.nanstd(Xdur[xcond]),
                    'pvel':np.nanmean(Xpvel[xcond]),'pvel_min':np.nanmin(Xpvel[xcond]),'pvel_max':np.nanmax(Xpvel[xcond]),'pvel_std':np.nanstd(Xpvel[xcond]),
                    'pvel2':np.nanmean(Xpvel2[xcond]),'pvel2_min':np.nanmin(Xpvel2[xcond]),'pvel2_max':np.nanmax(Xpvel2[xcond]),'pvel2_std':np.nanstd(Xpvel2[xcond])})

    ttt.name=streddy

    Table=Table.join(ttt,how='right')



    'Fiuguras'
    fsz=14;

    pc=plt.figure(figsize=(16.82,  3.63))

    gs = gridspec.GridSpec(1,11)
    ax1 = plt.subplot(gs[0, 0:2])

    if polx==-1:
        Xc,Yc=-Xamp*100,Xvel*100
        X=Xc[(Xc<=-vamp)&(Wc>=vrad)];Y=Yc[(Xc<=-vamp)&(Wc>=vrad)];
    else:
        Xc,Yc=Xamp*100,Xvel*100
        X=Xc[(Xc>=vamp)&(Wc>=vrad)];Y=Yc[(Xc>=vamp)&(Wc>=vrad)];


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
    plt.xlim()



    ax2 = plt.subplot(gs[0, 3:5])

    if polx==-1:
        Xc,Yc=-Xamp*100,Xeke
        X=Xc[(Xc<=-vamp)&(Wc>=vrad)];Y=Yc[(Xc<=-vamp)&(Wc>=vrad)];
    else:
        Xc,Yc=Xamp*100,Xeke
        X=Xc[(Xc>=vamp)&(Wc>=vrad)];Y=Yc[(Xc>=vamp)&(Wc>=vrad)];




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
    plt.xlim()


    ax3 = plt.subplot(gs[0, 6:8])

    if polx==-1:
        Xc,Yc=-Xamp*100,Xrad/1000
        X=Xc[(Xc<=-vamp)&(Wc>=vrad)];Y=Yc[(Xc<=-vamp)&(Wc>=vrad)];
    else:
        Xc,Yc=Xamp*100,Xrad/1000
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
    plt.ylim(40,Y.max())
    plt.xlim()


    ax4 = plt.subplot(gs[0,9:11])

    # X=Xc[Xc>=np.nanmean(Xc)-np.nanstd(Xc)];Y=Yc[Xc>=np.nanmean(Xc)-np.nanstd(Xc)];

    if polx==-1:
        xxc=-Xamp*100
        Xc,Yc=-Xrad/1000,Xpvel*100
        X=Xc[xxc<=-vamp];Y=Yc[xxc<=-vamp];
    else:
        xxc=Xamp*100
        Xc,Yc=Xrad/1000,Xpvel*100
        X=Xc[xxc>=vamp];Y=Yc[xxc>=vamp];



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
    plt.ylim(0,Y.max())
    plt.xlim()

    plt.tight_layout()
    plt.savefig('/home/igor/Documents/AntSimi_adt/figures/Scatter_{0}.png'.format(streddy))
    plt.close()



    'Figuras 3a,3b,3c - Frequencies (Formation, Duration Monthly /Yearly)'

    mostr=[u'Jan', u'Feb', u'Mar', u'Apr',u'May',u'Jun',u'Jul',u'Aug',u'Sep',u'Oct',u'Nov',u'Dec']


    plt.rcParams['font.weight']='bold'
    plt.rcParams['font.size']=12


    Bx=Xamp*100

    xmo = pd.DatetimeIndex(Xdate).month
    xfmo= pd.DatetimeIndex(Xfdate).month
    xyea= pd.DatetimeIndex(Xdate).year

    xmo=xmo[Bx>vamp];xyea=xyea[Bx>vamp];xfmo=xfmo[Bx>vamp]


    Xbox=np.zeros([12,len(np.arange(1993, dt.year.max() +1 ) )])

    for mi,yi in zip(xmo,xyea):
        Xbox[mi-1,yi-1993]=len(xmo[(xmo==mi)&(xyea==yi)])

    Xbox=pd.DataFrame(Xbox,columns=np.arange(1993, dt.year.max() +1 ),index=mostr)

    mrange=np.arange(0.5,13.5)
    yrange=np.arange(dt[0].year-0.5,dt[-1].year+1 +0.5)

    'BOX PLOT CYC ANTCYC'

    p=plt.figure(figsize=(13.59,  4.81 ))

    cbp=plt.boxplot(Xbox,widths=0.5,patch_artist=True)

    for box in cbp['boxes']:
        box.set(color='dimgray',linewidth=2)
        box.set(facecolor =colorx)
        # box.set(hatch = '/')

    plt.setp(cbp['whiskers'],color='dimgray',linewidth=2)
    plt.setp(cbp['fliers'],color='dimgray',linewidth=2)
    plt.setp(cbp['medians'],color='dimgray',linewidth=2)
    plt.setp(cbp['caps'],color='dimgray',linewidth=2)

    # for whisk in cbp
    #     whisk.set(color='dimgray',linewidth=2)
    #     # fly.set(color='k', linewidth=2,alpha=0.8)


    plt.grid()
    plt.xticks(ticks=np.arange(1,13),labels=mostr)
    # plt.ylim(5,35)

    plt.tight_layout()
    plt.savefig('/home/igor/Documents/AntSimi_adt/figures/Boxplot_{0}.png'.format(streddy))
    plt.close()

    'PDF annual frequency'
    p=plt.figure(figsize=(11.77,  3.86 ))

    n=plt.hist(xyea,bins=yrange,color=colorx,edgecolor='k',alpha=0.5,density=1)
    # plt.plot(np.arange(cyea[0],cyea[-1]+1),n[0],'-',color='navy',alpha=1)
    # sns.distplot(cyea,bins=yrange,color='navy',hist_kws={'edgecolor':'k'})
    ax1.set_xlabel('Years',fontweight='bold')

    ax1.set_xlim(dt[0].year-1,dt[-1].year+1)

    ax1.set_xticks(np.arange(dt[0].year,dt[-1].year+1)[::2])
    ax1.set_xticklabels(np.arange(dt[0].year,dt[-1].year+1)[::2])


    plt.tight_layout()
    plt.savefig('/home/igor/Documents/AntSimi_adt/figures/Histogram_frequency_anual_{0}.png'.format(streddy))
    plt.close()


Table.to_csv()

savep='/home/igor/Documents/AntSimi_adt/figures/Eddies_properties.dat'
Table.to_csv(savep,sep='\t',mode='w',float_format='%.6f',index=True,header=True,na_rep='NaN')



'''
Trajectories
'''


'CYCLONES'
ltmn1=-28;ltmx1=-12;lnmn1=-46;lnmx1=-34.

# lnmn=-43;lnmx=-36;ltmn=-25;ltmx=-18;
p=plt.figure(figsize=(8.33, 9.55))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmn1,llcrnrlat=ltmn1,urcrnrlat=ltmx1,urcrnrlon=lnmx1,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn1,lnmx1,2),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn1,ltmx1,2),labels=lbs,zorder=101,linewidth=0.1)




for itt in range(len(edatabase)):
    ltmn,ltmx,lnmn,lnmx,nnc,pole,colorx,streddy=edatabase[itt]


    if pole==-1:
        xtrack=nnc['track'][:]
        xlon,xlat=nnc['longitude'][:]-360,nnc['latitude'][:]

        ntrack=np.empty(xtrack.shape)*np.nan; #newtrack

        for ix in range(xtrack.max()):
            xts=xlon[xtrack==ix][0];
            yts=xlat[xtrack==ix][0];

            cond=(xts>=lnmn)&(xts<=lnmx)&(yts<=ltmx)&(yts>=ltmn)    

            if cond:
                ntrack[xtrack==ix]=xtrack[xtrack==ix]

        uni=np.unique(ntrack);uni=np.int_(uni[~np.isnan(uni)]);


        for vic in uni:

                xtt=nnc['longitude'][ntrack==vic]-360
                ytt=nnc['latitude'][ntrack==vic]

                xrr=np.mean(nnc['effective_radius'][:].data[ntrack==vic])/1e3
                xpp=np.mean(nnc['amplitude'][:].data[ntrack==vic])*100

                if (xrr>=40) & (xpp>=vamp):

                    # ncc=len(tracked[vic]['centers'])
                    # cp=tracked[vic]['polarity']
                    
                    plt.plot(xtt[0],ytt[0],color='k',marker='o',ms=1)
                    plt.plot(xtt,ytt,color=colorx,linewidth=0.5)
                    # plt.plot(*tracked[vic]['centers'][-1],color='k',marker='x',ms=3)

C.contourf(xxb,yyb,zb,colors='dimgray',alpha=0.5,levels=[-100,0],latlon=True)
C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=1,latlon=True)
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi_adt/figures/Trajectories_cyclones.png')
plt.close();



'ANTICYCLONES'
ltmn1=-28;ltmx1=-12;lnmn1=-46;lnmx1=-34.

# lnmn=-43;lnmx=-36;ltmn=-25;ltmx=-18;
p=plt.figure(figsize=(8.33, 9.55))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmn1,llcrnrlat=ltmn1,urcrnrlat=ltmx1,urcrnrlon=lnmx1,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn1,lnmx1,2),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn1,ltmx1,2),labels=lbs,zorder=101,linewidth=0.1)




for itt in range(len(edatabase)):
    ltmn,ltmx,lnmn,lnmx,nnc,pole,colorx,streddy=edatabase[itt]


    if pole==1:
        xtrack=nnc['track'][:]
        xlon,xlat=nnc['longitude'][:]-360,nnc['latitude'][:]

        ntrack=np.empty(xtrack.shape)*np.nan; #newtrack

        for ix in range(xtrack.max()):
            xts=xlon[xtrack==ix][0];
            yts=xlat[xtrack==ix][0];

            cond=(xts>=lnmn)&(xts<=lnmx)&(yts<=ltmx)&(yts>=ltmn)    

            if cond:
                ntrack[xtrack==ix]=xtrack[xtrack==ix]

        uni=np.unique(ntrack);uni=np.int_(uni[~np.isnan(uni)]);


        for vic in uni:

                xtt=nnc['longitude'][ntrack==vic]-360
                ytt=nnc['latitude'][ntrack==vic]

                xrr=np.mean(nnc['effective_radius'][:].data[ntrack==vic])/1e3
                xpp=np.mean(nnc['amplitude'][:].data[ntrack==vic])*100

                if (xrr>=40) & (xpp>=vamp):

                    # ncc=len(tracked[vic]['centers'])
                    # cp=tracked[vic]['polarity']
                    
                    plt.plot(xtt[0],ytt[0],color='k',marker='o',ms=1)
                    plt.plot(xtt,ytt,color=colorx,linewidth=0.5)
                    # plt.plot(*tracked[vic]['centers'][-1],color='k',marker='x',ms=3)

C.contourf(xxb,yyb,zb,colors='dimgray',alpha=0.5,levels=[-100,0],latlon=True)
C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=1,latlon=True)
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi_adt/figures/Trajectories_anticyclones.png')
plt.close();
