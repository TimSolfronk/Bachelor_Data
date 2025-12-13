import numpy as np

# model parameters in exaseis coordinate system
X_min = -13.5e3
X_max = 13.5e3
Y_min = 0
Y_max = 27.0e3
Z_min = -13.5e3
Z_max = 13.5e3

# generate input points for cvmh data query (in general)
dgrid = 500 # 500 m resolution
x = np.arange(X_min,X_max+dgrid,dgrid)
y = np.arange(Y_min,Y_max+dgrid,dgrid)
z = np.arange(Z_min,Z_max+dgrid,dgrid)

#print(x)
#print(y)
#print(z)

with open("cvmh_in_off-fault", "w") as text_file:
    for k in range (z.size):
        for j in range (y.size):
            for i in range (x.size):
                print('{} {} {}'.format(( 648446 - (0.5802386 * z[k]) - (0.8144465 * x[i])),
                                        (3625237 + (0.8144465 * z[k]) - (0.5802386 * x[i])), 
                                        np.maximum(y[j],100.0)), file=text_file)
                
fgrid = 100 # the fault has 100 m resolution
y_fault = np.arange(Y_min,Y_max+fgrid,fgrid)
z_fault = np.arange(Z_min,Z_max+fgrid,fgrid)
with open("cvmh_in_on-fault", "w") as text_file:
    for k in range (z_fault.size):
        for j in range (y_fault.size):
            print('{} {} {}'.format(( 648446 - (0.5802386 * z_fault[k]) - (0.8144465 * 0)),
                                    (3625237 + (0.8144465 * z_fault[k]) - (0.5802386 * 0)), 
                                    np.maximum(y_fault[j],100.0)), file=text_file)