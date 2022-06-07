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



bvar=Dataset('/home/igor/Documents/Extras/Sumida/MERCATOR/bat/etopo1_bedrock.nc')
lonb=bvar['lon'][:];latb=bvar['lat'][:];bat=bvar['Band1'][:]
lonb,latb=np.meshgrid(lonb,latb)


svar=Dataset('/home/igor/Documents/AntSimi/woa/woa18_decav_s00_04.nc')
tvar=Dataset('/home/igor/Documents/AntSimi/woa/woa18_decav_t00_04.nc')



T=tvar['t_an'][0,:,:,:]
S=svar['s_an'][0,:,:,:]

mask=T.mask

T=T.data;T[mask]=np.nan;
S=S.data;S[mask]=np.nan;


lon=tvar['lon'][:].data
lat=tvar['lat'][:].data
dep=tvar['depth'][:].data

lon,lat=np.meshgrid(lon,lat)

N=np.zeros(T.shape)*np.nan
H=np.zeros(T[0,:,:].shape)*np.nan

for yi in range(T.shape[1]):
    for xi in range(T.shape[2]):

        N[:-1,yi,xi]=np.sqrt(sw.bfrq(S[:,yi,xi],T[:,yi,xi],dep,lat=lat[yi,xi])[0].T)

        try:
            dind=np.argwhere(np.isnan(T[:,yi,xi])).T[0][0]-1 
            H[yi,xi]=dep[dind]

        except:
            H[yi,xi]=np.nan


H[mask[0,:,:]]=np.nan;




Rdi=(N*H/abs(sw.f(lat)))/1e3













