import scipy.signal as sg
import numpy as np
import matplotlib.pylab as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def smoo(self,nw):
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


def fig1(sizeh,sizev,lnmn,lnmx,ltmn,ltmx,grid_a=0.01):
    plt.rcParams['font.size']=18
    plt.rcParams['font.sans-serif']='Arial'
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.weight']='bold'


    fig = plt.figure(figsize=(sizeh,sizev))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent([lnmn, lnmx, ltmn, ltmx], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND,edgecolor='#fecf78')
    ax.add_feature(cfeature.LAKES)


    # ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.STATES)


    ax = plt.axes(projection=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True,alpha=grid_a)
    gl.top_labels = False
    gl.right_labels = False
    
def ecdf(data):
    n = len(data)
    x = np.sort(data)
    y = np.arange(1, n+1) / float(n)
    return x, y