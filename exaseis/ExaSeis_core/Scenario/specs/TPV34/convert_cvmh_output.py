import numpy as np
from netCDF4 import Dataset
import xarray as xr

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
outfile = open("cvmh_out",'r')
VpVsRho = []
for line in outfile:
    columns = (line.strip().split())
    # add vp, vs and rho, but already scale them correctly for exaseis
    VpVsRho.append([float(columns[Vp])/1000.0,
                          float(columns[Vs])/1000.0,
                          float(columns[Rho])/1000.0])
VpVsRho = np.asarray(VpVsRho)
# enforce minimum parameter values
VpVsRho[:,0][VpVsRho[:,0]<2.9840]=2.9840 # replace the Vp values less than 2984 m/s
VpVsRho[:,1][VpVsRho[:,1]<1.4000]=1.4000 # replace the Vs values less than 1400 m/s
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
rho = VpVsRho[:,:,:,2]                              # rho = rho
mu  = VpVsRho[:,:,:,1] ** 2 * rho                   # mu = Vs² * rho
lam = VpVsRho[:,:,:,0] ** 2 * rho - 2 * mu          # lam = Vp² * rho - 2 * mu

print('Writing NetCDF file')
nc = Dataset('tpv34_rhomulambda.nc', "w", format="NETCDF4")
nc.createDimension("x", len(x_nc))
nc.createDimension("y", len(y_nc))
nc.createDimension("z", len(z_nc))

vx = nc.createVariable("x","f4",("x",))
vx[:] = x_nc
vy = nc.createVariable("y","f4",("y",))
vy[:] = y_nc
vz = nc.createVariable("z","f4",("z",))
vz[:] = z_nc

mattype4 = np.dtype([('rho','f4'),('mu','f4'),('lambda','f4')])
mattype8 = np.dtype([('rho','f8'),('mu','f8'),('lambda','f8')])
mattype = nc.createCompoundType(mattype4,'material')

# transform to an array of tuples
arr = np.stack((rho,mu,lam), axis=3)
arr = arr.view(dtype=mattype8)              # format rho-mu-lam to tuples
#print(arr)
arr = arr.reshape(arr.shape[:-1])
print(arr[0,10,:])
print(arr[:,10,0])
mat = nc.createVariable("data",mattype,("x","y","z"))
mat[:] = arr
nc.close()
print('Done')
ds = xr.open_dataset('tpv34_rhomulambda.nc')
print(ds)