

import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from glob import glob


'Directory + Loading Data'


fnames = glob('/home/igor/Dropbox/META_eddydetection_images/data/*.pkl')
tnames = glob('/home/igor/Dropbox/META_eddydetection_images/T/*.pkl')
xnames = glob('/home/igor/Dropbox/META_eddydetection_images/xy/X_*.pkl')
ynames = glob('/home/igor/Dropbox/META_eddydetection_images/xy/Y_*.pkl')

fnames.sort()
tnames.sort()
xnames.sort()
ynames.sort()



data = {}
for fname in fnames:
    eddy = fname.split('/')[-1].split('.')[0]
    D = pd.read_pickle(fname).values

    D = xr.DataArray(D,dims=['eddy','day'],
                     coords={
                         'eddy':np.arange(D.shape[0]),
                         'day':np.arange(D.shape[1])
                     })
    data[eddy] = D

tdata = {}
for tname in tnames:
    teddy = tname.split('/')[-1].split('.')[0]
    T = pd.read_pickle(tname).values

    T = xr.DataArray(T,dims=['time','day'],
                     coords={
                         'time':np.arange(T.shape[0]),
                         'day':np.arange(T.shape[1])
                     })
    tdata[teddy] = T

xdata = {}
for xname in xnames:
    xeddy = xname.split('/')[-1].split('.')[0]
    X = pd.read_pickle(xname).values

    X = xr.DataArray(X,dims=['x','pts','day'],
                     coords={
                         'x':np.arange(X.shape[0]),
                         'pts':np.arange(X.shape[1]),
                         'day':np.arange(X.shape[2])
                     })
    xdata[xeddy] = X

ydata = {}
for yname in ynames:
    yeddy = yname.split('/')[-1].split('.')[0]
    Y = pd.read_pickle(yname).values

    Y = xr.DataArray(Y,dims=['x','pts','day'],
                     coords={
                         'x':np.arange(Y.shape[0]),
                         'pts':np.arange(Y.shape[1]),
                         'day':np.arange(Y.shape[2])
                     })
    ydata[yeddy] = Y




'Funcion setting'
avg = lambda D,std: D.mean('eddy')+std*D.std('eddy')
gwr = lambda D,std: (avg(D,std).differentiate('day')/avg(D,std))
std = lambda D: np.abs(np.vstack([gwr(D,-1)-gwr(D,0),gwr(D,1)-gwr(D,0)])).max(0)

func = lambda x,a,b,c: a*np.exp(-x/b)


def get_fit(D):
    fit = []
    for gwri in (D.differentiate('day')/D).values:
        popt,_=curve_fit(func, np.arange(len(gwri))[~np.isnan(gwri)], gwri[~np.isnan(gwri)])
        fit.append(func(np.arange(len(gwri)),*popt))
    fit = np.vstack(fit)*xr.ones_like(D)
    fit = fit.where((fit.min('day')>=0)&(fit.max('day')<=0.4),drop=True)

    return fit





'Plotting Growth Rate'


plt.rcParams['font.size']=20
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.weight']='bold'


kw = {
    'xlim':[0,20],
    'ylim':[-0.015,0.2]
}
    
for E in data:


    ne=2

    fig,ax = plt.subplots()

    fit = get_fit(data[E])
    x = data[E].day
    y = (data[E].differentiate('day')/data[E]).mean('eddy')
    popt,_=curve_fit(func, x[~np.isnan(y)],y[~np.isnan(y)])
    ax.plot(data[E].day,func(data[E].day,*popt),lw=2.5,color='0.2',zorder=1)
    fitstd = fit.std('eddy')
    ax.fill_between(fit.day,func(data[E].day,*popt)-fitstd,func(data[E].day,*popt)+fitstd,color='0.9',zorder=0)
    (data[E].differentiate('day')/data[E]).mean('eddy').plot(ax=ax,marker='o',color='#78ca9cff',lw=0,markeredgecolor='0.3',zorder=5)
    ax.axvline(ne*popt[1],color='0.2',linestyle='--')
    ax.text(ne*popt[1]+1,0.16,'{:.1f}$\,$days'.format(ne*popt[1]),color='0.3')
    ax.text(3,-0.01,'N = {}'.format(data[E].shape[0]),ha='center')
    
    ax.set(
        title=E,
        xlim=[0,20],
        ylim=[-0.1,0.20],
    )
    print(E)
    print(popt[1])
    print( func(data[E].day,*popt)[:5].data)
    


'Plotting Mean Growth Rate'
colorcode=['crimson','mediumseagreen','cornflowerblue','lightsalmon','chocolate','darkorange','navy']
ixd=0
for E in data:
    fig,ax = plt.subplots()
#     ax.scatter(data[E].day*xr.ones_like(data[E]),data[E],0.1,c='0.2')
#     _ = data[E].plot.line(x='day',color='0.2',add_legend=False,lw=0.1)
    _ = data[E].rolling(day=5,min_periods=1).mean().plot.line(x='day',color='0.2',add_legend=False,lw=0.1)
    _ = data[E].mean('eddy').plot.line(x='day',color=colorcode[ixd],add_legend=False,lw=2)
    ax.set(
        title=E,
        xlim=[1,30],
        ylim=[20,120],
    )
    ixd+=1
    plt.grid()
    plt.xlabel('Time [day]')
    plt.ylabel('Radius [km]')

    plt.tight_layout()
    plt.savefig('/home/igor/Documents/AntSimi/figures/fullgrowth_{0}.png'.format(E))
    plt.close()


'Plotting Mean + All growth rates'

fig,ax = plt.subplots(figsize=(11,8.45))
ixd=0
for E in data:
    _ = data[E].mean('eddy').plot.line(x='day',alpha=0.5,add_legend=False,lw=2,color=colorcode[ixd],label=E)

    ixd+=1

ax.legend(ncol=4)
ax.set(
    xlim=[1,20],
    ylim=[30,75],
)
plt.ylim(30,90)
plt.xlim(0,25)
plt.grid()
plt.xlabel('Time [day]')
plt.ylabel('Radius [km]')

plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/Mean_growth_rate.png')
plt.close()


'Radius x Time'

E = 'RCE'
# plt.close('all')
# for i in range(100):
i=36
fig,ax = plt.subplots(figsize=(6.62, 4.32))
data[E].isel(eddy=i).rolling(day=2,min_periods=1).mean().plot(ax=ax,color='chocolate',linewidth=3)
ax.set(title=E,xlim=[0,50],ylim=[20,100])
plt.grid()
plt.ylabel(r'Radius [km]')
plt.xlabel(r'Time [day]')

plt.xlim(0,25)
plt.ylim(0,120)
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/Isolated_growth_rate_RCE_{0}.png'.format(i))
plt.close()





'Radius x Time'

E = 'AE'
# plt.close('all')
# for i in range(100):
i=52
fig,ax = plt.subplots(figsize=(6.62, 4.32))
data[E].isel(eddy=i).rolling(day=2,min_periods=1).mean().plot(ax=ax,color='crimson',linewidth=3)
ax.set(title=E,xlim=[0,50],ylim=[20,100])
plt.grid()
plt.ylabel(r'Radius [km]')
plt.xlabel(r'Time [day]')

plt.xlim(0,25)
plt.ylim(0,120)
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/Isolated_growth_rate_{0}_{1}.png'.format(E,i))
plt.close()




E = 'CFE'
# plt.close('all')
# for i in range(100):
i=7
fig,ax = plt.subplots(figsize=(6.62, 4.32))
data[E].isel(eddy=i).rolling(day=2,min_periods=2).mean().plot(ax=ax,color='mediumseagreen',linewidth=3)
ax.set(title=E,xlim=[0,50],ylim=[20,100])
plt.grid()
plt.ylabel(r'Radius [km]')
plt.xlabel(r'Time [day]')

plt.xlim(1,25)
plt.ylim(0,120)
plt.tight_layout()

plt.savefig('/home/igor/Documents/AntSimi/figures/Isolated_growth_rate_CSTE.png')
plt.close()


'Radius Growth - Horizontal'

from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset

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

'Topography'
ltmn=-29;ltmx=-9;lnmn=-52;lnmx=-32.

topo=Dataset('/home/igor/Documents/Extras_Works/Sumida/HYCOM/bat/etopo1_bedrock.nc')
xb=topo['lon'][:];yb=topo['lat'][:];zb=topo['Band1'][:];
zb=smoo(zb,20)
xb,yb,zb=geocut(xb,yb,zb,ltmn,lnmn,ltmx,lnmx);
xxb,yyb=np.meshgrid(xb,yb);


'Royal Charlotte Eddy'

'limits'
ltmn1=-20;ltmx1=-14;lnmn1=-41;lnmx1=-36.


for i in range(xdata['X_'+E].shape[0]):
    E = 'RCE'

    p=plt.figure(figsize=(6.6 , 7.23))
    lbs=[True, False,False,True]
    C = Basemap(llcrnrlon=lnmn1,llcrnrlat=ltmn1,urcrnrlat=ltmx1,urcrnrlon=lnmx1,resolution='h',projection='cyl')

    C.drawcoastlines(linewidth=1);
    C.fillcontinents(color='#fecf78')
    C.drawstates(linewidth=0.5)
    C.drawmeridians(np.arange(lnmn1+1,lnmx1,2),labels=lbs,zorder=100,linewidth=0.2)
    C.drawparallels(np.arange(ltmn1,ltmx1,2),labels=lbs,zorder=101,linewidth=0.2)

    # i=30
    # alphax=0.1

    # for dx in range(100):
    #     plt.plot(xdata['X_'+E][i,:,dx],ydata['Y_'+E][i,:,dx],color='chocolate',alpha=alphax)

    #     alphax+=0.3
    zz=(np.zeros(ydata['Y_'+E][i,:,:25].shape)+1)*np.arange(1,26)[None,:]
    plt.contour(xdata['X_'+E][i,:,:25],ydata['Y_'+E][i,:,:25],zz,cmap='inferno_r',levels=np.arange(1,26),zorder=27)

    C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,0],latlon=True,zorder=25)
    C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=1,latlon=True,zorder=26)
    plt.tight_layout()
    plt.savefig('/home/igor/Documents/AntSimi/figures/growth_maps/horizontal_growth_{0}_{1}.png'.format(E,i))
    plt.close()



'Cabo Frio Eddy'

'limits'
ltmn1=-28;ltmx1=-20;lnmn1=-46;lnmx1=-39.

# lnmn=-43;lnmx=-36;ltmn=-25;ltmx=-18;



for i in range(xdata['X_'+E].shape[0]):

    E = 'CFE'
    p=plt.figure(figsize=(8.83, 7.23))
    lbs=[True, False,False,True]
    C = Basemap(llcrnrlon=lnmn1,llcrnrlat=ltmn1,urcrnrlat=ltmx1,urcrnrlon=lnmx1,resolution='h',projection='cyl')

    C.drawcoastlines(linewidth=1);
    C.fillcontinents(color='#fecf78')
    C.drawstates(linewidth=0.5)
    C.drawmeridians(np.arange(lnmn1+1,lnmx1,2),labels=lbs,zorder=100,linewidth=0.1)
    C.drawparallels(np.arange(ltmn1,ltmx1,2),labels=lbs,zorder=101,linewidth=0.1)
    # i=7

    # colorx=0.01

    # alphax=0.05
    # for dx in range(25):
    #     plt.plot(xdata['X_'+E][i,:,dx],ydata['Y_'+E][i,:,dx],color='mediumseagreen',alpha=alphax)
    #     alphax+=0.03
    zz=(np.zeros(ydata['Y_'+E][i,:,:25].shape)+1)*np.arange(1,26)[None,:]
    plt.contour(xdata['X_'+E][i,:,:25],ydata['Y_'+E][i,:,:25],zz,cmap='inferno_r',levels=np.arange(1,26))
    # plt.colorbar() 

    C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,0],latlon=True,zorder=25)
    C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=1,latlon=True,zorder=26)
    plt.tight_layout()
    plt.savefig('/home/igor/Documents/AntSimi/figures/growth_maps/horizontal_growth_{0}_{1}.png'.format(E,i))
    plt.close()






'Propagation Speed Case Study'
'Map'

E = 'RCE'
i=29
ltmn1=-20;ltmx1=-14;lnmn1=-41;lnmx1=-36.


import seawater as sw
days=25
xx=xdata['X_'+E][i,:,:days].mean('pts')
yy=ydata['Y_'+E][i,:,:days].mean('pts')


'Modo 1'
dist=sw.dist(lon=xx,lat=yy)[0]
k=np.insert(dist,0,0)
dist1=np.cumsum(dist)
pvel=((dist1[-1]*1e3)/(days*86400))*100

'Modo 1'
dist2=sw.dist(lon=[xx[0],xx[-1]],lat=[yy[0],yy[-1]])[0]
pvel2=((dist2[-1]*1e3)/(days*86400))*100

'Modo 3'
dist=sw.dist(lon=xx.mean(),lat=yy)[0]
k=np.insert(dist,0,0)
dist3=np.cumsum(dist)[-1]
pvel3=((dist3*1e3)/(days*86400))*100

'Modo 4'
dist4=sw.dist(lon=xx.mean(),lat=[yy[0],yy[-1]])[0][0]
pvel4=((dist4*1e3)/(days*86400))*100

'Modo 5'
dist5=sw.dist(lon=xx,lat=yy)[0].mean()
pvel5=((dist5*1e3)/(86400))*100





p=plt.figure(figsize=(8.83, 7.23))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmn1,llcrnrlat=ltmn1,urcrnrlat=ltmx1,urcrnrlon=lnmx1,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn1+1,lnmx1,2),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn1,ltmx1,2),labels=lbs,zorder=101,linewidth=0.1)
# i=7

zz=(np.zeros(ydata['Y_'+E][i,:,:25].shape)+1)*np.arange(1,26)[None,:]
plt.contour(xdata['X_'+E][i,:,:25],ydata['Y_'+E][i,:,:25],zz,cmap='inferno_r',levels=np.arange(1,26),zorder=27)

plt.text(-36.5,-14.,r'd1={0:.1f},c1= {1:.2f}'.format(dist1[-1],pvel),ha='center')
plt.text(-36.5,-14.4,r'd2={0:.1f},c2= {1:.2f}'.format(dist2[0],pvel2),ha='center')
plt.text(-36.5,-14.8,r'd3={0:.1f},c3= {1:.2f}'.format(dist3,pvel3),ha='center')
plt.text(-36.5,-15.2,r'd4={0:.1f},c4= {1:.2f}'.format(dist4,pvel4),ha='center')
plt.text(-36.5,-15.4,r'c5= {0:.2f}'.format(pvel5),ha='center')



# plt.scatter(xdata['X_'+E][i,:,:25].mean('pts'),ydata['Y_'+E][i,:,:25].mean('pts'),c=zz[0,:],cmap='inferno_r')
# plt.colorbar() 

C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,0],latlon=True,zorder=25)
C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=1,latlon=True,zorder=26)
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/Example_Propagation_{0}_{1}.png'.format(E,i))
plt.close()





E = 'CFE'
i=7
ltmn1=-28;ltmx1=-20;lnmn1=-46;lnmx1=-39.


import seawater as sw
days=25
xx=xdata['X_'+E][i,:,:days].mean('pts')
yy=ydata['Y_'+E][i,:,:days].mean('pts')
xx[xx<-50]=np.nan 
yy[yy>-5]=np.nan

'Modo 1'
dist=sw.dist(lon=xx,lat=yy)[0]
k=np.insert(dist,0,0)
dist1=np.nancumsum(dist)
pvel=((dist1[-1]*1e3)/(days*86400))*100

'Modo 1'
dist2=sw.dist(lon=[xx[0],xx[-1]],lat=[yy[0],yy[-1]])[0]
pvel2=((dist2[-1]*1e3)/(days*86400))*100

'Modo 3'
dist=sw.dist(lon=xx.mean(),lat=yy)[0]
k=np.insert(dist,0,0)
dist3=np.nancumsum(dist)[-1]
pvel3=((dist3*1e3)/(days*86400))*100


'Modo 5'
dist5=np.nanmean(sw.dist(lon=xx,lat=yy)[0])
pvel5=((dist5*1e3)/(86400))*100





p=plt.figure(figsize=(8.83, 7.23))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmn1,llcrnrlat=ltmn1,urcrnrlat=ltmx1,urcrnrlon=lnmx1,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn1+1,lnmx1,2),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn1,ltmx1,2),labels=lbs,zorder=101,linewidth=0.1)
# i=7

zz=(np.zeros(ydata['Y_'+E][i,:,:25].shape)+1)*np.arange(1,26)[None,:]
plt.contour(xdata['X_'+E][i,:,:25],ydata['Y_'+E][i,:,:25],zz,cmap='inferno_r',levels=np.arange(1,26),zorder=27)

plt.text(-44.5,-26,r'd1={0:.1f},c1= {1:.2f}'.format(dist1[-1],pvel),ha='center')
plt.text(-44.5,-26.4,r'd2={0:.1f},c2= {1:.2f}'.format(dist2[0],pvel2),ha='center')
plt.text(-44.5,-26.8,r'd3={0:.1f},c3= {1:.2f}'.format(dist3,pvel3),ha='center')
plt.text(-44.5,-27.2,r'd4={0:.1f},c4= {1:.2f}'.format(dist4,pvel4),ha='center')
plt.text(-44.5,-27.6,r'c5= {0:.2f}'.format(pvel5),ha='center')



# plt.scatter(xdata['X_'+E][i,:,:25].mean('pts'),ydata['Y_'+E][i,:,:25].mean('pts'),c=zz[0,:],cmap='inferno_r')
# plt.colorbar() 

C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,0],latlon=True,zorder=25)
C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=1,latlon=True,zorder=26)
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/Example_Propagation_{0}_{1}.png'.format(E,i))
plt.close()







E = 'AE'
i=52
ltmn1=-20;ltmx1=-14;lnmn1=-41;lnmx1=-36.


import seawater as sw
days=11
xx=xdata['X_'+E][i,:,:days].mean('pts')
yy=ydata['Y_'+E][i,:,:days].mean('pts')
xx[xx<-50]=np.nan 
yy[yy>-5]=np.nan


'Modo 1'
dist=sw.dist(lon=xx,lat=yy)[0]
k=np.insert(dist,0,0)
dist1=np.nancumsum(dist)
pvel=((dist1[-1]*1e3)/(days*86400))*100

'Modo 1'
dist2=sw.dist(lon=[xx[0],xx[-1]],lat=[yy[0],yy[-1]])[0]
pvel2=((dist2[-1]*1e3)/(days*86400))*100

'Modo 3'
dist=sw.dist(lon=xx.mean(),lat=yy)[0]
k=np.insert(dist,0,0)
dist3=np.nancumsum(dist)[-1]
pvel3=((dist3*1e3)/(days*86400))*100

'Modo 4'
dist4=sw.dist(lon=xx.mean(),lat=[yy[0],yy[-1]])[0][0]
pvel4=((dist4*1e3)/(days*86400))*100

'Modo 5'
dist5=sw.dist(lon=xx,lat=yy)[0].mean()
pvel5=((dist5*1e3)/(86400))*100





p=plt.figure(figsize=(8.83, 7.23))
lbs=[True, False,False,True]
C = Basemap(llcrnrlon=lnmn1,llcrnrlat=ltmn1,urcrnrlat=ltmx1,urcrnrlon=lnmx1,resolution='h',projection='cyl')

C.drawcoastlines(linewidth=1);
C.fillcontinents(color='#fecf78')
C.drawstates(linewidth=0.5)
C.drawmeridians(np.arange(lnmn1+1,lnmx1,2),labels=lbs,zorder=100,linewidth=0.1)
C.drawparallels(np.arange(ltmn1,ltmx1,2),labels=lbs,zorder=101,linewidth=0.1)
# i=7

zz=(np.zeros(ydata['Y_'+E][i,:,:days].shape)+1)*np.arange(1,days+1)[None,:]
# plt.contour(xdata['X_'+E][i,:,:days][0],ydata['Y_'+E][i,:,:days][0],zz[0],cmap='inferno_r',levels=np.arange(1,days+1),zorder=27)

xx1=xdata['X_'+E][i,:,:days]
yy1=ydata['Y_'+E][i,:,:days]
xx1=np.array(xx1)
yy1=np.array(yy1)
xx1[xx1<-50]=np.nan 
yy1[yy1>-5]=np.nan

plt.plot(xx1[:,0],yy1[:,0],'orange',zorder=27,label=r't=0')
plt.plot(xx1[:,-1],yy1[:,-1],color='purple',zorder=27,label=r't=11')

plt.legend(loc='upper left')
plt.text(-36.5,-14.,r'd1={0:.1f},c1= {1:.2f}'.format(dist1[-1],pvel),ha='center')
plt.text(-36.5,-14.4,r'd2={0:.1f},c2= {1:.2f}'.format(dist2[0],pvel2),ha='center')
plt.text(-36.5,-14.8,r'd3={0:.1f},c3= {1:.2f}'.format(dist3,pvel3),ha='center')
plt.text(-36.5,-15.2,r'd4={0:.1f},c4= {1:.2f}'.format(dist4,pvel4),ha='center')
# plt.text(-36.5,-15.4,r'c5= {0:.2f}'.format(pvel5),ha='center')



# plt.scatter(xdata['X_'+E][i,:,:25].mean('pts'),ydata['Y_'+E][i,:,:25].mean('pts'),c=zz[0,:],cmap='inferno_r')
# plt.colorbar() 

C.contourf(xxb,yyb,zb,colors='dimgray',alpha=1,levels=[-100,0],latlon=True,zorder=25)
C.contour(xxb,yyb,-zb,colors='k',alpha=1,levels=[0,100],linewidths=1,latlon=True,zorder=26)
plt.tight_layout()
plt.savefig('/home/igor/Documents/AntSimi/figures/Example_Propagation_{0}_{1}.png'.format(E,i))
plt.close()