import numpy as np

expansion_size = 1
zero = 1.e-4
fractional_coord = True # If the basis is given in units of the primitive vectors (fractional coordinates)

primitive_vectors = [ #Angstrom
#  [ 5.6649999619,         0.0000000000,         0.0000000000],
#  [-2.8324999809,         4.9060338794,         0.0000000000],
#  [ 0.0000000000,         0.0000000000,         4.5310001373],
 [1.0,  0.0             ,     0.0], 
 [0.0,  1.73205080756888,     0.0], 
 [0.0,  0.0             ,     1.0], 
]

base_vectors = np.array([ # Fractional coordinates
#  [ 0.838800013,         0.161199987,        0.250000000],
#  [ 0.322400004,         0.161200002,        0.250000000],
#  [ 0.838800013,         0.677600026,        0.250000000],
#  [ 0.161199987,         0.322399974,        0.750000000],
#  [ 0.677600026,         0.838800013,        0.750000000],
#  [ 0.161199987,         0.838800013,        0.750000000],
 [-0.382000,   0.3820000,   0.750000],
 [ 0.382000,   0.3820000,   0.250000],
 [-0.118000,   0.1180000,   0.250000],
 [ 0.118000,   0.1180000,   0.750000],
 [-0.118000,  -0.1180000,   0.250000],
 [ 0.118000,  -0.1180000,   0.750000],
 [-0.382000,  -0.3820000,   0.750000],
 [ 0.382000,  -0.3820000,   0.250000],
 [ 0.0000  ,   0.33333  ,   0.50000],
 [ 0.0000  ,  -0.33333  ,   0.50000], 
 [ 0.0000  ,   0.33333  ,   0.00000],
 [ 0.0000  ,  -0.33333  ,   0.00000], 
 [-0.5000  ,   0.16667  ,   0.50000],
 [-0.5000  ,  -0.16667  ,   0.50000],
 [-0.5000  ,   0.16667  ,   0.00000],
 [-0.5000  ,  -0.16667  ,   0.00000],
])

Jij_list =[
# distances were obtained using Vesta and the file: Mn3Sn/extra/0.199_Mn3Sn.vasp
# [Rij (Ang), Jij (meV), Dijz (meV), different unit cell],
[ 2.76302, -3.61, 0.00], # J1
[ 2.92541,  5.36, 0.52], # J2
[ 2.73959,-11.77, 0.19], # J3
[ 4.02396,  8.62, 0.00], # J4
[ 3.89097,  9.40, 0.00], # J5
[ 4.53100,  1.47, 0.00], # J6
[ 4.92003, -4.92, 0.00], # J7
[ 4.90691, -6.43, 0.00], # J8
# [ 1.5811388300841898, -3.89, 0.00], # J9
# [ 1.5811388300841898, -7.99, 0.00], # J10
# [ 1.5811388300841898,  0.92, 0.00], # J11
]

a1 = np.array(primitive_vectors[0])
a2 = np.array(primitive_vectors[1])
a3 = np.array(primitive_vectors[2])


nn_vect_list = [
#    [base_vectors[1] - base_vectors[0],  1.0],                        # W  arrow J2
#    [base_vectors[2] - base_vectors[0], -1.0],                        # NW arrow J2
#    [base_vectors[1] - base_vectors[0] + np.array([1, 0, 0]), -1.0],  # E  arrow J3 
#    [base_vectors[2] - base_vectors[0] + np.array([0,-1, 0]),  1.0]   # SE arrow J3
]

rotation_list = [
    0* 2.0*np.pi/6.0,
    4* 2.0*np.pi/6.0,
    2* 2.0*np.pi/6.0,
    5* 2.0*np.pi/6.0,
    1* 2.0*np.pi/6.0, 
    3* 2.0*np.pi/6.0,
]

from scipy.spatial.transform import Rotation as R
rotation_axis = np.array([0, 0, 1])

# vec = [1,0,0]
# rotation_angle = 6 * 2.0*np.pi/6.0
# rotation = R.from_rotvec(rotation_angle * rotation_axis)
# rotated_vec = rotation.apply(vec)
# print(rotated_vec)


print("{:3}{:3} {:3}{:3}{:3}  {:14}{:14}{:14}{:14}   ".format("  i", "  j", "  da", " db", " dc", "   Dijx", "   Dijy", "   Dijz", "   Jij") )

# Finding pair for each atom in the basis
for atomi, R_i in enumerate(base_vectors) :
    # if atomi in [0,1,2,3,4,5] :
        if fractional_coord :
            R_i = R_i[0]*a1 +  R_i[1]*a2 + R_i[2]*a3
        
        # Rotation to correct identify the DMI orientation
        rotation_angle = rotation_list[atomi]
        rotation = R.from_rotvec(rotation_angle * rotation_axis)
        # print('angle degree', rotation_angle/np.pi * 180)

        # Loop over atoms in the origin and surrounding unit cells
        for i in range(-expansion_size, expansion_size+1):
            for j in range(-expansion_size, expansion_size+1):
                for k in range(-expansion_size, expansion_size+1):
                    for atomj, R_j in enumerate(base_vectors) :
                        if fractional_coord :
                            R_j = R_j[0]*a1 +  R_j[1]*a2 + R_j[2]*a3

                        R_j = i*a1 + j*a2 + k*a3 + R_j
                        distance = np.linalg.norm( np.array(R_i) - np.array(R_j) )

                        for [Jij_distance, Jij, Dij] in Jij_list :
                            if abs(distance - Jij_distance) < zero :
                                
                                Dij_aux = 0.0
                                if Dij > zero :
                                    
                                    for [R_k, chirality] in nn_vect_list :
                                        if fractional_coord :
                                            R_k = R_k[0]*a1 +  R_k[1]*a2 + R_k[2]*a3
                                                                     
                                        rotated_vec = rotation.apply(R_k)
                                        # print(rotated_vec, np.array(R_j) - np.array(R_i))

                                        if np.linalg.norm( np.array(R_j) - np.array(R_i) - rotated_vec  ) < zero :
                                            # print("DMI")
                                            Dij_aux = chirality*Dij
                                            break
                                    

                                print("{:3}{:3} {:3}{:3}{:3}  {:14.8f}{:14.8f}{:14.8f}{:14.8f}  {:.16f}".format(atomi, atomj, i, j, k, 0., 0., Dij_aux, Jij, distance) )

                                  