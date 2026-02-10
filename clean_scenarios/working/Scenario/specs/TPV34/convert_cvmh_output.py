import numpy as np
from netCDF4 import Dataset
import xarray as xr
import sys

# domain size parameters
X_min = -16.2
X_max = 16.2
Y_min = 0
Y_max = 32.4
Z_min = -16.2
Z_max = 16.2

# read last 3 columns from the cvmh output files
Vp = 17-1; Vs = 18-1; Rho = 19-1 

### for the general domain (just get material properties)
outfile = open("cvmh_out_off-fault",'r')
rhoVsVp = []
for line in outfile:
    columns = (line.strip().split())
    # add vp, vs and rho, but already scale them correctly for exaseis
    rhoVsVp.append([float(columns[Rho])/1000.0,
                    float(columns[Vs])/1000.0,
                    float(columns[Vp])/1000.0])
    
rhoVsVp = np.asarray(rhoVsVp)
# enforce minimum parameter values
rhoVsVp[:,0][rhoVsVp[:,0]<2.22034]=2.22034 # replace the Rho values less than 2220.34 m/s
rhoVsVp[:,1][rhoVsVp[:,1]<1.4]=1.4 # replace the Vs values less than 1400 m/s
rhoVsVp[:,2][rhoVsVp[:,2]<2.984]=2.984 # replace the Vp values less than 2984 m/s


### write NetCDF files
# generate tpv34_rhovsvp.nc file
dgrid = 0.6
x_nc = np.arange(X_min,X_max+dgrid,dgrid)
y_nc = np.arange(Y_min,Y_max+dgrid,dgrid)
z_nc = np.arange(Z_min,Z_max+dgrid,dgrid)

rhoVsVp = rhoVsVp.reshape((len(z_nc),len(y_nc),len(x_nc),3))

print('Writing NetCDF file for whole domain')
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

mattype4 = np.dtype([('rho','f4'),('vs','f4'),('vp','f4')])
mattype8 = np.dtype([('rho','f8'),('vs','f8'),('vp','f8')])
mattype = nc.createCompoundType(mattype4,'material')

arr = rhoVsVp
arr = arr.view(dtype=mattype8)
arr = arr.reshape(arr.shape[:-1])

mat = nc.createVariable("data",mattype,("z","y","x"))

mat[:] = arr
nc.close()
print('Done')



### for the fault (get mu/mu0 multiplier)
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


VsRho = VsRho.reshape((len(y_fault),len(z_fault),2))
mu0 = 32038120320.0 # in Pa

mu_mult = np.zeros((y_fault.size, z_fault.size, 3))

# compute mu/mu0 and add the y and z coordinate as extra elements to access them more easily in the yaml file 
# (saving y and z here does not lead to inaccuracies, as the values are interpolated linearly between entries in the yaml file)
for j in range (y_fault.size):
    for k in range (z_fault.size):
        mu_mult[j,k]  = [y_fault[j],z_fault[k],VsRho[j,k,0] ** 2 * VsRho[j,k,1] / mu0]


print('Writing NetCDF file for fault')
nc = Dataset('tpv34_mu_mult_plus.nc', "w", format="NETCDF4")
nc.createDimension("y", len(y_fault))
nc.createDimension("z", len(z_fault))

vy = nc.createVariable("y","f4",("y",))
vy[:] = y_fault
vz = nc.createVariable("z","f4",("z",))
vz[:] = z_fault

mattype4 = np.dtype([('y_coord','f4'),('z_coord','f4'),('mu_mult','f4')])
mattype8 = np.dtype([('y_coord','f8'),('z_coord','f8'),('mu_mult','f8')])
mattype = nc.createCompoundType(mattype4,'material')

arr = mu_mult

arr = arr.view(dtype=mattype8)           
mat = nc.createVariable("data",mattype,("y","z"))
mat[:] = arr
nc.close()
print('Done')
