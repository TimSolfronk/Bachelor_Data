import numpy as np
from netCDF4 import Dataset

# model parameters
X_min = -13.5e3
X_max = 13.5e3
Y_min = 0
Y_max = 27.0e3
Z_min = -13.5e3
Z_max = 13.5e3

fault_width_min = 0
fault_width_max = 15e3
fault_length_min = -15e3
fault_length_max = 15e3

# read last 3 columns from the 'outfile' files
Vp = 17-1; Vs = 18-1; Rho = 19-1 # python begin from 0!

