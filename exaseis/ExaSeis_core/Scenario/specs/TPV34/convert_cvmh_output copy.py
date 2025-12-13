import numpy as np
from netCDF4 import Dataset
import xarray as xr
import sys

# model parameters
X_min = -13.5
X_max = 13.5
Y_min = 0
Y_max = 27.0
Z_min = -13.5
Z_max = 13.5

fault_width_min = 0
fault_width_max = 15.0

# read last 3 columns from the 'outfile' files
Vp = 17-1; Vs = 18-1; Rho = 19-1 # python begin from 0!

### inside refinement zone
outfile = open("cvmh_out_off-fault",'r')
VpVsRho = []
for line in outfile:
    columns = (line.strip().split())
    # add vp, vs and rho, but already scale them correctly for exaseis
    VpVsRho.append([float(columns[Vp])/1000.0,
                          float(columns[Vs])/1000.0,
                          float(columns[Rho])/1000.0])
VpVsRho = np.asarray(VpVsRho)
# enforce minimum parameter values
VpVsRho[:,0][VpVsRho[:,0]<2.984]=2.984 # replace the Vp values less than 2984 m/s
VpVsRho[:,1][VpVsRho[:,1]<1.4]=1.4 # replace the Vs values less than 1400 m/s
VpVsRho[:,2][VpVsRho[:,2]<2.22034]=2.22034 # replace the Rho values less than 2220.34 m/s


### write NetCDF files
# generate tpv34_rhomulambda-inner.nc file
dgrid = 0.5
x_nc = np.arange(X_min,X_max+dgrid,dgrid)
y_nc = np.arange(Y_min,Y_max+dgrid,dgrid)
z_nc = np.arange(Z_min,Z_max+dgrid,dgrid)
#print(x_nc)
#print(y_nc)
#print(z_nc)


VpVsRho = VpVsRho.reshape((len(z_nc),len(y_nc),len(x_nc),3))

print('Writing NetCDF file')
nc = Dataset('tpv34_rhovsvp.nc', "w", format="NETCDF4")
nc.createDimension("x", len(x_nc))
nc.createDimension("y", len(y_nc))
nc.createDimension("z", len(z_nc))

vz = nc.createVariable("z","f4",("z",))
vz[:] = z_nc
vy = nc.createVariable("y","f4",("y",))
vy[:] = y_nc
vx = nc.createVariable("x","f4",("x",))
vx[:] = x_nc

mattype4 = np.dtype([('vs','f4'),('vp','f4'),('rho','f4')])
mattype8 = np.dtype([('vs','f8'),('vp','f8'),('rho','f8')])
mattype = nc.createCompoundType(mattype4,'material')

# transform to an array with 4 dim (x,y,z,0-2)
arr = VpVsRho
arr = arr.view(dtype=mattype8)              # format to dtype
#print(arr)
arr = arr.reshape(arr.shape[:-1])
print(arr[0,0,0])
mat = nc.createVariable("data",mattype,("z","y","x"))
mat[:] = arr
nc.close()
print('Done')







outfile = open("cvmh_out_on-fault",'r')
VsRho = []
for line in outfile:
    columns = (line.strip().split())
    VsRho.append([float(columns[Vs]),
                          float(columns[Rho])])
VsRho = np.asarray(VsRho)
# enforce minimum parameter values
VsRho[:,0][VsRho[:,0]<1400.0]=1400.0 # replace the Vs values less than 1400 m/s
VsRho[:,1][VsRho[:,1]<2220.34]=2220.34 # replace the Rho values less than 2220.34 m/s


fgrid = 0.1 # the fault has 100 m resolution
y_fault = np.arange(Y_min,Y_max+fgrid,fgrid)
z_fault = np.arange(Z_min,Z_max+fgrid,fgrid)


VsRho = VsRho.reshape((len(z_fault),len(y_fault),2))
mu0 = 32038120320.0 # in Pa
mu_mult  = VsRho[:,:,0] ** 2 * VsRho[:,:,1] / mu0       # mu_mult = mu/mu0   (with mu = Vs² * rho)     (shear modulus)

print('Writing NetCDF file')
nc = Dataset('tpv34_mu_mult.nc', "w", format="NETCDF4")
nc.createDimension("z", len(z_fault))
nc.createDimension("y", len(y_fault))

vy = nc.createVariable("y","f4",("y",))
vy[:] = y_fault
vz = nc.createVariable("z","f4",("z",))
vz[:] = z_fault

mattype4 = np.dtype([('mu_mult','f4')])
mattype8 = np.dtype([('mu_mult','f8')])
mattype = nc.createCompoundType(mattype4,'material')

# transform to an array with 4 dim (y,z,0-1)
arr = mu_mult

arr = arr.view(dtype=mattype8)              # format rho-mu to dtype
#print(arr)
print(arr[0,0])
mat = nc.createVariable("data",mattype,("z","y"))
mat[:] = arr
nc.close()
print('Done')

y_coord =  10.0
z_coord =  10.0 

with xr.open_dataset('tpv34_rhovsvp.nc') as ds:
    print(ds)
    print()
    print()
    print(str(round(ds.sel(x=0.0,y=y_coord,z=z_coord,method="nearest")["data"].item()[2],5)) + ", " + str(round(ds.sel(x=0.0,y=y_coord,z=z_coord,method="nearest")["data"].item()[1] ** 2 * ds.sel(x=0.0,y=y_coord,z=z_coord,method="nearest")["data"].item()[2]* 1e9/mu0,8)))
with xr.open_dataset('tpv34_mu_mult.nc') as ds2:
    print(ds2)
    print(ds2.sel(y=y_coord,z=z_coord, method="nearest")["data"].values)