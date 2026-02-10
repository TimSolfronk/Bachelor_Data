import numpy as np
from netCDF4 import Dataset

# fault size
Y_min = 0
Y_max = 15.5
Z_min = -20.0
Z_max = 20.0

# read last 2 columns from tpv35_input_data.txt
mu_s = 5-1; s_xz = 6-1; 

### inside refinement zone
outfile = open("tpv35_input_data.txt",'r')
muS_sXZ = []
for line in outfile:
    columns = (line.strip().split())
    muS_sXZ.append([float(columns[mu_s]),
                          float(columns[s_xz])])
muS_sXZ = np.asarray(muS_sXZ)


### write NetCDF files
dgrid = 0.1
y_nc = np.arange(Y_min,Y_max+dgrid/2.0,dgrid)
z_nc = np.arange(Z_min,Z_max+dgrid/2.0,dgrid)


muS_sXZ = muS_sXZ.reshape((len(y_nc),len(z_nc),2))

print('Writing NetCDF file')
nc = Dataset('tpv35_muS_sXZ.nc', "w", format="NETCDF4")
nc.createDimension("y", len(y_nc))
nc.createDimension("z", len(z_nc))

vy = nc.createVariable("y","f4",("y",))
vy[:] = y_nc
vz = nc.createVariable("z","f4",("z",))
vz[:] = z_nc

mattype4 = np.dtype([('mu_s','f4'),('s_xz','f4')])
mattype8 = np.dtype([('mu_s','f8'),('s_xz','f8')])
mattype = nc.createCompoundType(mattype4,'material')

arr = muS_sXZ
arr = arr.view(dtype=mattype8)  # format to tuples on lowest level instead of another array dim
arr = arr.reshape(arr.shape[:-1])
mat = nc.createVariable("data",mattype,("y","z"))
mat[:] = arr
nc.close()
print('Done')

'''
# sample reading of netcdf file

import xarray as xr

y_coord = 0.0
z_coord = -20.0
with xr.open_dataset('tpv35_muS_sXZ.nc') as ds:
    print(ds)
    print()
    print(ds.sel(y=y_coord,z=z_coord, method="nearest")["data"].values)
'''
