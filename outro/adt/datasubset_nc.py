import xarray as xr
from glob import glob

path='/home/igor/Documents/AntSimi/altdata/'
files=glob(path+'*.nc')
files.sort()

for ix in range(len(files)):
    data = xr.open_dataset(files[ix])

    for i in range(data.time.shape[0]):
        day = data.isel(time=i)
        day.to_netcdf('/home/igor/Documents/AntSimi_adt/data/'+str(day.time.data).split('T')[0].replace('-','_')+'.nc', format='NETCDF4')
