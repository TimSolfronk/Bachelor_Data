import numpy as np
from netCDF4 import Dataset
import xarray as xr

# model parameters
Y_min = 0
Y_max = 25.0
Z_min = -25.0
Z_max = 25.0

fault_width_min = 0
fault_width_max = 15.0

# read last 3 columns from the 'outfile' files
y = 4-1; x = 5-1; # python begin from 0!

### inside refinement zone
outfile = open("tpv29_tpv30_geometry_ext_25m_data.txt",'r')
x_data = []
for line in outfile:
    columns = (line.strip().split())
    if(len(columns) < 7):
        continue
    if(float(columns[y])<0):
        continue
    # add the data, but already scale them correctly for exaseis
    x_data.append(float(columns[x])/1000.0)
x_data = np.asarray(x_data)

#print(input_data[-3:-1,:])
### write NetCDF files
# generate tpv16_fault.nc file
dgrid = 0.025
y_nc = np.arange(Y_min,Y_max+dgrid/2.0,dgrid)
z_nc = np.arange(Z_min,Z_max+dgrid/2.0,dgrid)


x_data = x_data.reshape((len(y_nc),len(z_nc),1))
#print(input_data[:,0,:])

print('Writing NetCDF file')
nc = Dataset('tpv29_fault_disp.nc', "w", format="NETCDF4")
nc.createDimension("y", len(y_nc))
nc.createDimension("z", len(z_nc))

vy = nc.createVariable("y","f4",("y",))
vy[:] = y_nc
vz = nc.createVariable("z","f4",("z",))
vz[:] = z_nc

mattype4 = np.dtype([('x','f4')])
mattype8 = np.dtype([('x','f8')])
mattype = nc.createCompoundType(mattype4,'material')

arr = x_data
arr = arr.view(dtype=mattype8)  # format to tuples on lowest level instead of another array dim
arr = arr.reshape(arr.shape[:-1])
mat = nc.createVariable("data",mattype,("y","z"))
mat[:] = arr
nc.close()
print('Done')

tracers = [
    [15.6,  2.0],
    [ 5.0,  5.0],
    [12.0,  5.0],
    [ 1.4,  9.0],
    [10.1, 11.1],
    [ 0.0, 15.0], #surface above fault hypocenter
    [10.0, 15.0], #on fault hypocenter
    [16.0, 15.0],
    [ 6.1, 15.8],
    [12.0, 20.0],
    [ 6.2, 15.7],
    [ 6.0, 15.4],
    [ 0.0, 25.0],
    [12.0, 25.0],
    [ 5.7, 25.1],
    [ 4.7, 25.9],
    [15.3, 29.0],
    [15.3, 29.3],
    [ 5.0, 30.0],
    [11.0, 30.0],
    [ 0.0, 35.0],
    [13.0, 35.0],
    [10.5, 36.7],
    [ 4.5, 37.0],
]
y_coord = 15.6
z_coord = 2.0
with xr.open_dataset('tpv29_buffer_n.nc') as ds2:
    print(ds2)
    print(ds2.sel(y=y_coord * 1000,x=z_coord * 1000, method="nearest")["z"].values)
with xr.open_dataset('tpv29_fault_disp.nc') as ds2:
    print(ds2)
    print(ds2.sel(y=y_coord,z=z_coord, method="nearest")["data"].values)
    print()
    for coords in tracers:
       value = float(str(ds2.sel(y=coords[0],z=coords[1]-20.0, method="nearest")["data"].values)[1:-2])
       print("[" + str(round(19.8 + value,3)) + ", "+ str(coords[0]) + ", " + str(coords[1]) + "], \t# x_offset = " + str(value)) 
    print()
    for coords in tracers:
       value = float(str(ds2.sel(y=coords[0],z=coords[1]-20.0, method="nearest")["data"].values)[1:-2])
       print("[" + str(round(20.2 + value,3)) + ", "+ str(coords[0]) + ", " + str(coords[1]) + "], \t# x_offset = " + str(value)) 
