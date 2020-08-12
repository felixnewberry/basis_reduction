from fenics import *
from mshr import *
import numpy as np
#import dolfin as dfn
import scipy.io as sciio
import scipy.io
import scipy.sparse as sparse
from scipy.io import savemat
from scipy.io import loadmat
#import os
import time

#set_log_level(WARNING)
#set_log_active(False)

################################################################################
### LDC Approach

print("run_LDC_ensemble underway")
#Record simulation time
start = time.time()

from NavierStokesLDC import Navier_Stokes_LDC

# Load Random Variables

#content = loadmat('./LDC_data/u_nu_vec_2.mat')
content = loadmat('./LDC_data/u_nu_vec_hi.mat')
# content = loadmat('./LDC_data/u_nu_vec.mat')

u_lid_vec = content['u_top_vec']
nu_vec = content['nu_vec']



# content = loadmat('./sample_i.mat')
# sample_i = int(content['sample_i'])

# content = loadmat('./LDC_data/inputs_vec.mat')
content = loadmat('./LDC_data/inputs_vec_high.mat')

nx = int(content['nx'])
delta_u = float(content['delta_u'])
delta_nu = float(content['delta_nu'])
# delta_nu_1 = float(content['delta_nu_1'])

#run_count = content['run_count']
#QoI = float(content['QoI'])

## TO obtain high-fidelity data:
# nx = int(64)

# Note on timing: in practice - to compare model expense the solve is timed, not the loading to high-fidelity mesh and interpolation - this could be mvoed to a post-processing step
# Time 200 low: 2.49s
# Time 10 high 19.32s
# Time ratio 155.58

#n_samps = 200
n_samps = 2000

u_matrix_0 = np.zeros((n_samps,65))
u_matrix_1 = np.zeros((n_samps,65))
u_matrix_2 = np.zeros((n_samps,65))
u_matrix_3 = np.zeros((n_samps,65))
u_matrix_4 = np.zeros((n_samps,65))

# matrix
u_matrix_5 = np.zeros((n_samps,4225))
u_matrix_6 = np.zeros((n_samps,4225))
u_matrix_7 = np.zeros((n_samps,4225))


# Apply deltas to u and nu % want them to be multiplicative. :/
nu_vec = nu_vec*(1+delta_nu)
u_lid_vec = u_lid_vec*(1+delta_u)

#print(QoI)



for i in range(n_samps): #200

    sample_i = i;
    #sample_i = int(Nan_vec[0,i])
    u_lid = float(u_lid_vec[sample_i])
    u_lid = Constant(u_lid)
    #
    nu = float(nu_vec[sample_i])
    #nu = nu*(1+delta_nu)
    #nu_1 = nu*(1+delta_nu_1)
    nu = Constant(nu)
    #
    #u_y_array, u_x_array, p_array_mid, p_array_vert, p_array_base = Navier_Stokes_LDC(u_lid, nu_0, nu_1, nx)

    # # linearly varying nu and sigmoid
    u_y_array, u_x_array, p_array_mid, p_array_vert, p_array_base, p_field, u_x_field, u_y_field = Navier_Stokes_LDC(u_lid, nu, nx)
    #
    u_matrix_0[i,:] = u_y_array
    u_matrix_1[i,:] = u_x_array
    u_matrix_2[i,:] = p_array_mid
    u_matrix_3[i,:] = p_array_vert
    u_matrix_4[i,:] = p_array_base

    u_matrix_5[i,:] = p_field
    u_matrix_6[i,:] = u_x_field
    u_matrix_7[i,:] = u_y_field
    aa = Navier_Stokes_LDC(u_lid, nu, nx)

# Save QoI

# scipy.io.savemat('./u_meshes/u_matrix.mat', mdict={'u_matrix_0':u_matrix_0, 'u_matrix_1':u_matrix_1, 'u_matrix_2':u_matrix_2, 'u_matrix_3':u_matrix_3, 'u_matrix_4':u_matrix_4, 'u_matrix_5':u_matrix_5, 'u_matrix_6':u_matrix_6, 'u_matrix_7':u_matrix_7})

scipy.io.savemat('u_matrix_hi.mat', mdict={'u_matrix_0':u_matrix_0, 'u_matrix_1':u_matrix_1, 'u_matrix_2':u_matrix_2, 'u_matrix_3':u_matrix_3, 'u_matrix_4':u_matrix_4, 'u_matrix_5':u_matrix_5, 'u_matrix_6':u_matrix_6, 'u_matrix_7':u_matrix_7})

end = time.time()
print(end-start)
