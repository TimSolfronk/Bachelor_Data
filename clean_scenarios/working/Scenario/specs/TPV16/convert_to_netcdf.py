import numpy as np
from netCDF4 import Dataset
import xarray as xr

# model parameters
Y_min = 0
Y_max = 19.5
Z_min = 0.0
Z_max = 48.0

fault_width_min = 0
fault_width_max = 15.0

# read last all data columns of 'tpv16_input_file.txt'
# data is saved at:
sig_n = 5-1; sig_h_shear = 6-1; sig_v_shear = 7-1; mu_s = 10-1; mu_d = 11-1; d_c = 12-1; cohesion = 13-1; time_forced_rup = 14-1

outfile = open("tpv16_input_file.txt",'r')
input_data = []
for line in outfile:
    columns = (line.strip().split())
    if(len(columns) < 14):
        continue
    # add the data, but already scale them correctly for exaseis
    input_data.append([-float(columns[sig_n])/1000.0/1000.0,
                          float(columns[sig_h_shear])/1000.0/1000.0,
                          float(columns[sig_v_shear])/1000.0/1000.0,
                          float(columns[mu_s]),
                          float(columns[mu_d]),
                          float(columns[d_c]),
                          float(columns[cohesion])/1000.0/1000.0,
                          float(columns[time_forced_rup])])
input_data = np.asarray(input_data)


### write NetCDF files
# generate tpv16_fault.nc file  (75m resolution)
dgrid = 0.075
y_nc = np.arange(Y_min,Y_max+dgrid/2.0,dgrid)
z_nc = np.arange(Z_min,Z_max+dgrid/2.0,dgrid)


input_data = input_data.reshape((len(y_nc),len(z_nc),8))

print('Writing NetCDF file')
nc = Dataset('tpv16_fault.nc', "w", format="NETCDF4")
nc.createDimension("y", len(y_nc))
nc.createDimension("z", len(z_nc))

vy = nc.createVariable("y","f4",("y",))
vy[:] = y_nc
vz = nc.createVariable("z","f4",("z",))
vz[:] = z_nc

mattype4 = np.dtype([('s_xx','f4'),('s_xz','f4'),('s_xy','f4'),('mu_s','f4'),('mu_d','f4'),('d_c','f4'),('cohesion','f4'),('forced_rupture_time','f4')])
mattype8 = np.dtype([('s_xx','f8'),('s_xz','f8'),('s_xy','f8'),('mu_s','f8'),('mu_d','f8'),('d_c','f8'),('cohesion','f8'),('forced_rupture_time','f8')])
mattype = nc.createCompoundType(mattype4,'material')

arr = input_data
arr = arr.view(dtype=mattype8)  # format to tuples on lowest level instead of another array dim     
arr = arr.reshape(arr.shape[:-1])      
mat = nc.createVariable("data",mattype,("y","z"))
mat[:] = arr
nc.close()
print('Done')



'''
# sample reading of netcdf file

import xarray as xr

y_coord = 19.5
z_coord = 48.0
with xr.open_dataset('tpv16_fault.nc') as ds:
    print(ds)
    print()
    print(ds.sel(y=y_coord,z=z_coord, method="nearest")["data"].values)
'''
