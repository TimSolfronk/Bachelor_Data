import numpy as np

# model parameters
X_min = -13.5e3
X_max = 13.5e3
Y_min = 0
Y_max = 27.0e3
Z_min = -13.5e3
Z_max = 13.5e3


fault_width_min = 0
fault_width_max = 15e3

# generate input points for cvmh data query (refinement zone)
dgrid = 500 # the inner zone has 500 m resolution
x = np.arange(X_min/2,X_max/2+dgrid,dgrid)
y = np.arange(Y_min/2,Y_max/2+dgrid,dgrid)
z = np.arange(fault_width_min,fault_width_max+dgrid,dgrid)
with open("cvmh_in", "w") as text_file:
    for k in range (z.size):
        for j in range (y.size):
            for i in range (x.size):
                print('{} {} {}'.format(( 648446 - (0.5802386 * x[i]) - (0.8144465 * y[j])),
                                        (3625237 + (0.8144465 * x[i]) - (0.5802386 * y[j])), 
                                        np.maximum(z[k],100.0)), file=text_file)
                
