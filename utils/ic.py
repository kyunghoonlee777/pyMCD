import numpy as np

def arccos(x):
    return np.arccos(x)

def sin(x):
    return np.sin(x)

def cos(x):
    return np.cos(x)

def get_vector(xyz_coordinates,idx1,idx2,normalize = False):
    start = xyz_coordinates[idx1]
    end = xyz_coordinates[idx2]
    v = end-start
    if normalize:
        d = np.linalg.norm(v)
        v /= d
    return v

def get_cross_vector(v1,v2,normalize = False):
    v = np.cross(v1,v2)
    if normalize:
        d = np.linalg.norm(v)
        v /= d
    return v

def get_distance(xyz_coordinates,idx1,idx2):
    v = get_vector(xyz_coordinates,idx1,idx2)
    return np.linalg.norm(v)

def get_vector_angle(v1,v2):
    return np.arccos(np.dot(v1,v2))

def get_angle(xyz_coordinates,idx1,idx2,idx3):
    # Formula: theta = cos-1(u.v/|u||v|)
    v_21 = get_vector(xyz_coordinates,idx2,idx1,True)
    v_23 = get_vector(xyz_coordinates,idx2,idx3,True)
    return get_vector_angle(v_21,v_23)

def get_dihedral_angle(xyz_coordinates,idx1,idx2,idx3,idx4):
    # Formula: theta = cos-1((uxw)_n.(vxw)_n], (uxw)_n: normalized cross vector
    v_12 = get_vector(xyz_coordinates,idx1,idx2,True) # u
    v_23 = get_vector(xyz_coordinates,idx2,idx3,True) # w
    v_34 = get_vector(xyz_coordinates,idx3,idx4,True) # v
    w1 = get_cross_vector(v_12,v_23,True)
    w2 = get_cross_vector(v_23,v_34,True)
    return get_vector_angle(w1,w2)

def get_single_q_element(xyz_coordinates,internal_coordinate):
    if len(internal_coordinate) == 2:
        idx1,idx2 = internal_coordinate
        return get_distance(xyz_coordinates,idx2,idx1)
    elif len(internal_coordinate) == 3:
        idx1,idx2,idx3 = internal_coordinate
        return get_angle(xyz_coordinates,idx1,idx2,idx3)
    else:
        idx1,idx2,idx3,idx4 = internal_coordinate
        return get_dihedral_angle(xyz_coordinates,idx1,idx2,idx3,idx4)

def get_q(xyz_coordinates,internal_coordinates):
    q = []
    for internal_coordinate in internal_coordinates:
        q.append(get_single_q_element(xyz_coordinates,internal_coordinate))
    return np.array(q).T

# REF: Molecular Vibrations (1955), the theory of infrared and raman vibrational spectra, Wilson Bright, J. C. Decius, Paul C. cross
def get_wilsonB_vector(xyz_coordinates,internal_coordinate):
    vector_dict = dict()
    if len(internal_coordinate) == 2:
        idx1,idx2 = internal_coordinate
        e_12 = get_vector(xyz_coordinates,idx1,idx2,True)
        vector_dict[idx2] = e_12
        vector_dict[idx1] = -e_12
    elif len(internal_coordinate) == 3:
        idx1,idx2,idx3 = internal_coordinate
        e_21 = get_vector(xyz_coordinates,idx2,idx1,True)
        e_23 = get_vector(xyz_coordinates,idx2,idx3,True)
        r_12 = get_distance(xyz_coordinates,idx2,idx1)
        r_32 = get_distance(xyz_coordinates,idx3,idx2)
        theta = get_vector_angle(e_21,e_23)
        if np.abs(theta-np.pi) > 0.001:
            vector_dict[idx1] = (cos(theta)*e_21 - e_23)/(r_12*sin(theta))
            vector_dict[idx2] = ((r_12-r_32*cos(theta))*e_21 + (r_32-r_12*cos(theta))*e_23)/(r_12*r_32*sin(theta))
            vector_dict[idx3] = (cos(theta)*e_23 - e_21)/(r_32*sin(theta))
        else:
            vector_dict[idx1] = (cos(theta)*e_21 - e_23)/(r_12*sin(theta))
            vector_dict[idx3] = (cos(theta)*e_23 - e_21)/(r_32*sin(theta))
    else:
        idx1,idx2,idx3,idx4 = internal_coordinate
        r_21 = get_distance(xyz_coordinates,idx1,idx2)
        r_32 = get_distance(xyz_coordinates,idx3,idx2)
        r_43 = get_distance(xyz_coordinates,idx4,idx3)
        e_12 = get_vector(xyz_coordinates,idx1,idx2,True)
        e_23 = get_vector(xyz_coordinates,idx2,idx3,True)
        e_32 = get_vector(xyz_coordinates,idx3,idx2,True)
        e_43 = get_vector(xyz_coordinates,idx4,idx3,True)
        theta_2 = get_vector_angle(-e_12,e_23)
        theta_3 = get_vector_angle(e_32,-e_43)
        u = get_cross_vector(e_12,e_23,True)
        v = get_cross_vector(e_43,e_32,True)
        vector_dict[idx1] = -u/(r_21*sin(theta_2))
        vector_dict[idx2] = (r_32 - r_21*cos(theta_2))/(r_32*r_21*sin(theta_2))*u - cos(theta_3)/(r_32*sin(theta_3)) * v
        vector_dict[idx3] = (r_32 - r_43*cos(theta_3))/(r_32*r_43*sin(theta_3))*v - cos(theta_2)/(r_32*sin(theta_2)) * u
        vector_dict[idx4] = -v/(r_43*sin(theta_3))
    return vector_dict

def get_wilsonB_matrix(xyz_coordinates,internal_coordinates):
    # Reshape xyzs
    n = len(xyz_coordinates)
    m = len(internal_coordinates)
    B_matrix = np.zeros((m,3*n))
    for row_idx,internal_coordinate in enumerate(internal_coordinates):
        vectors = get_wilsonB_vector(xyz_coordinates,internal_coordinate)
        for index in vectors:
            B_matrix[row_idx,3*index:3*index+3] = vectors[index]
    return B_matrix

def get_general_inverse(B):
    G = B.T@B # matrix multiplication
    D,V = np.linalg.eigh(G)
    D = np.where(D>1e-4,1/D,0)
    D_inverse = np.diag(D) # General D_inverse
    G_inverse = V@D_inverse@V.T
    return G_inverse@B.T   

def solve_equation(B,delta_q):
    # Solve: delta_q = B*delta_x -> delta_x = (B.T*B)^(-1)*B.T*delta_q
    B_inverse = get_general_inverse(B)
    delta_x = B_inverse@delta_q
    return delta_x


def get_xyz_updates(xyz_coordinates,internal_coordinates,delta_q): # Note that here, delta_qs[internal_coordinate] = delta_q
    B = get_wilsonB_matrix(xyz_coordinates,internal_coordinates)
    delta_x = solve_equation(B,delta_q) # 3N x 1 -> N x 3
    n = len(xyz_coordinates)
    delta_x = np.reshape(delta_x,(n,3))
    return delta_x


def update_geometry(xyz_coordinates,q_updates,criteria = 1e-4,max_iteration = 30):
    internal_coordinates = list(q_updates.keys())
    delta_q = np.array(list(q_updates.values())).T
    q = get_q(xyz_coordinates,internal_coordinates)
    #print ('q: ',q)
    #print ('internals: ',internal_coordinates)
    target_q = q + delta_q
    original_coordinates = np.copy(xyz_coordinates)
    for i in range(max_iteration):
        delta_x = get_xyz_updates(xyz_coordinates,internal_coordinates,delta_q)
        xyz_coordinates += delta_x
        q = get_q(xyz_coordinates,internal_coordinates)
        delta_q = target_q - q
        x_norm = np.linalg.norm(delta_x)
        if x_norm < criteria:
            #print (delta_q,q_updates,q)
            break
        elif i == max_iteration-1:
            print (f'Iteration did not converged! The final updated vector has norm: {x_norm}')
            for xyz in original_coordinates:
                print (xyz)
            print ('update: ',q_updates)

def get_internal_coordinate_info(coordinate_list,internal_coordinates):
    infos = dict()
    q = get_q(coordinate_list,internal_coordinates)
    for idx,internal_coordinate in enumerate(internal_coordinates):
        infos[internal_coordinate] = q[idx]
    return infos


if __name__ == '__main__':
    from pyMCD import chem
    molecule = chem.Molecule('R.xyz')
    coordinate_list = molecule.get_coordinate_list()
    internal_coordinates = molecule.get_redundant_coordinates()
    #print (internal_coordinates)
    #exit()
    delta_qs = dict()
    delta_qs[(0,1)] = 0.0
    delta_qs[(1,2)] = 0.0
    delta_qs[(2,3)] = 0.0
    #delta_qs[(0,1,2,3)] = np.pi/3
    #delta_qs[(1,2,3)] = 0.0
    #delta_qs[(0,1,2,3)] = np.pi*2/3
    internal_coordinates = list(delta_qs.keys())
    print (get_q(coordinate_list,internal_coordinates))
    update_geometry(coordinate_list,delta_qs)
    for coordinate in coordinate_list:
        print (coordinate)
    print (get_internal_coordinate_info(coordinate_list,internal_coordinates))
