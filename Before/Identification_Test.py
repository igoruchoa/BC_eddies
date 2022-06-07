from netCDF4 import Dataset
from py_eddy_tracker.dataset.grid import RegularGridDataset, UnRegularGridDataset
# import logging
# logging.basicConfig(level=logging.DEBUG)
from datetime import *
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from glob import glob

'Rodando o dado de Eddy Identification, obs: melhor usar no pront'



'Station'
path='/home/igor/Documents/AntSimi/data/'
files=glob(path + '*.nc');files.sort()

for ix in range(len(files)):


    nc=Dataset(files[ix])
    lonc,latc=nc['longitude'][:],nc['latitude'][:];
    slac=nc['sla'][:]
    adtc=nc['adt'][:]

    dt0=datetime(1950,1,1)
    dt=dt0 + timedelta(days=np.int(nc['time'][:]))
    # ddt=pd.DatetimeIndex([dt0 + timedelta(days=np.int(dia)) for dia in nc['time'][:]])



    h = RegularGridDataset(files[ix], 'longitude', 'latitude')

    'Filter'
    h.high_filter('sla', 10, 5)
    # h.bessel_high_filter('sla', 500, order=3)

    ac, cc = h.eddy_identification('sla', 'ugosa', 'vgosa', date=dt,step= 0.0025)
    # ac, cc = h.eddy_identification('sla', 'ugosa', 'vgosa', date=ddt,step= 0.0025)



    with Dataset('/home/igor/Documents/AntSimi/tmp/Anticyclonic_{0}.nc'.format(dt.strftime('%Y%m%d') ), 'w') as h:
        ac.to_netcdf(h)
    with Dataset('/home/igor/Documents/AntSimi/tmp/Cyclonic_{0}.nc'.format(dt.strftime('%Y%m%d') ), 'w') as h:
        cc.to_netcdf(h)

    # with Dataset('/home/igor/Documents/AntSimi/tmp_bessel/Anticyclonic_{0}.nc'.format(dt.strftime('%Y%m%d') ), 'w') as h:
    #     ac.to_netcdf(h)
    # with Dataset('/home/igor/Documents/AntSimi/tmp_bessel/Cyclonic_{0}.nc'.format(dt.strftime('%Y%m%d') ), 'w') as h:
    #     cc.to_netcdf(h)




    fig = plt.figure(figsize=(9.65, 8.75))
    ax = fig.add_axes([.03,.03,.94,.94])
    ax.set_title('Eddies detected -- Anticyclonic(red) and Cyclonic(blue)')
    ax.set_ylim(latc.min(),latc.max())
    ax.set_xlim(lonc.min(),lonc.max())
    ax.set_aspect('equal')
    h.display(ax,'sla',cmap='seismic')
    # ax.contourf(lonc,latc,slac,cmap='seismic',levels=np.linspace(-0.5,0.5,100))
    # ax.contourf(lonc,latc,adtc[ix],cmap='jet',levels=np.linspace(0,1,100))

    ac.display(ax, color='m', linewidth=.5)
    cc.display(ax, color='m', linewidth=.5)
    # ax.grid()

    plt.tight_layout()
    fig.savefig('/home/igor/Downloads/testes_EddyT/AntSimi/Eddy_Id_{0}.png'.format(ix))
    plt.close()