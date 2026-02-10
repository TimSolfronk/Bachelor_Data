import numpy as np

# model parameters in exaseis coordinate system
X_min = -16.2e3
X_max = 16.2e3
Y_min = 0
Y_max = 32.4e3
Z_min = -16.2e3
Z_max = 16.2e3

# generate input points for cvmh data query (in general)
dgrid = 600 # resolution of 1/2 cell
x = np.arange(X_min,X_max+dgrid,dgrid)
y = np.arange(Y_min,Y_max+dgrid,dgrid)
z = np.arange(Z_min,Z_max+dgrid,dgrid)

with open("cvmh_in_off-fault", "w") as text_file:
    for k in range (z.size):
        for j in range (y.size):
            for i in range (x.size): # convert to utm coordinate system
                print('{} {} {}'.format(( 648446 - (0.5802386 * z[k]) - (0.8144465 * x[i])),
                                        (3625237 + (0.8144465 * z[k]) - (0.5802386 * x[i])), 
                                        np.maximum(y[j],100.0)), file=text_file)
                
fgrid = 100 # the fault has 100 m resolution (1/12 cell)
y_fault = np.arange(Y_min,Y_max+fgrid,fgrid)
z_fault = np.arange(Z_min,Z_max+fgrid,fgrid)
with open("cvmh_in_on-fault", "w") as text_file:
    for j in range (y_fault.size):
        for k in range (z_fault.size): # convert to utm coordinate system
            print('{} {} {}'.format(( 648446 - (0.5802386 * z_fault[k]) - (0.8144465 * 0)),
                                    (3625237 + (0.8144465 * z_fault[k]) - (0.5802386 * 0)), 
                                    np.maximum(y_fault[j],100.0)), file=text_file)
