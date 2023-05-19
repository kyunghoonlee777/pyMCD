import numpy as np
import itertools


# REF: Molecular Vibrations (1955), the theory of infrared and raman vibrational spectra, Wilson Bright, J. C. Decius, Paul C. cross
# REF: The efficient optimization of molecular geometries using redundant internal coordinates (2002), 

def arccos(x):
    return np.arccos(x)

def sin(x):
    return np.sin(x)

def cos(x):
    return np.cos(x)

def get_vector(xyz_coordinates,idx1,idx2,normalize = False):
    start = xyz_coordinates[idx1]
    end = xyz_coordinates[idx2]
    v = start-end
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

def kronecker_delta(a,b):
    if a==b:
        return 1
    else:
        return 0

def sign_factor(a,b,c):
    return kronecker_delta(a,b) - kronecker_delta(a,c)

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

def get_wilsonB_vector(xyz_coordinates,internal_coordinate):
    vector_dict = dict()
    if len(internal_coordinate) == 2:
        idx1,idx2 = internal_coordinate
        u = get_vector(xyz_coordinates,idx1,idx2,True)
        # a: idx1, m: idx1, n: idx2
        vector_dict[idx2] = sign_factor(idx2,idx1,idx2)*u
        vector_dict[idx1] = sign_factor(idx1,idx1,idx2)*u
    elif len(internal_coordinate) == 3:
        idx1,idx2,idx3 = internal_coordinate
        u = get_vector(xyz_coordinates,idx1,idx2,False)
        l_u = np.linalg.norm(u)
        u /= l_u # Normalized u
        v = get_vector(xyz_coordinates,idx3,idx2,False)
        l_v = np.linalg.norm(v) # Normalized v
        v /= l_v
        w = get_cross_vector(u,v,True) # Normalized w
        uw = get_cross_vector(u,w,False)
        wv = get_cross_vector(w,v,False)
        # m: idx1, o: idx2, n: idx3
        for idx in [idx1,idx2,idx3]:
            vector_dict[idx] = sign_factor(idx,idx1,idx2) * uw/l_u + sign_factor(idx,idx3,idx2) * wv/l_v
    else:
        idx1,idx2,idx3,idx4 = internal_coordinate
        # m:idx1, o:idx2, p:idx3, n: idx4
        u = get_vector(xyz_coordinates,idx1,idx2,False)
        l_u = np.linalg.norm(u)
        u /= l_u # Normalized u
        v = get_vector(xyz_coordinates,idx4,idx3,False)
        l_v = np.linalg.norm(v)
        v /= l_v # Normalized v
        w = get_vector(xyz_coordinates,idx3,idx2,False)
        l_w = np.linalg.norm(w)
        w /= l_w # Normalized w
        uw = get_cross_vector(u,w,False)
        vw = get_cross_vector(v,w,False)
        cos_phi_u = np.sum(u*w)
        sin_phi_u = np.sqrt(1-cos_phi_u**2)
        cos_phi_v = np.sum(v*w)
        sin_phi_v = np.sqrt(1-cos_phi_v**2)
        cos_q = np.sum(uw*vw)/(sin_phi_u*sin_phi_v)
        for idx in [idx1,idx2,idx3,idx4]:
            first_term = sign_factor(idx,idx1,idx2) * uw/(l_u*sin_phi_u**2)
            second_term = sign_factor(idx,idx3,idx4) * vw/(l_v*sin_phi_v**2)
            third_term = cos_phi_u*uw/(l_w*sin_phi_u**2) - cos_phi_v*vw/(l_w*sin_phi_v**2)
            third_term *= sign_factor(idx,idx2,idx3)
            vector_dict[idx] = first_term + second_term + third_term

    return vector_dict

def get_wilsonB_derivative_sub_matrix(xyz_coordinates,internal_coordinate):
    m = len(internal_coordinate)
    indices = list(range(m))
    scanning_indices = list(itertools.product(indices,indices))
    tensor = np.zeros((m,m,3,3))
    if m == 2:
        a,b = internal_coordinate
        u = get_vector(xyz_coordinates,a,b,False)
        l_u = np.linalg.norm(u)
        u /= l_u # Normalized u
        # m: a, n: b
        for indices in list(scanning_indices):
            x,y = indices
            idx1 = internal_coordinate[x]
            idx2 = internal_coordinate[y]
            # a: idx1, b: idx2
            for i in range(3):
                for j in range(3):
                    u_i = u[i]
                    u_j = u[j]
                    tensor[x][y][i][j] = (-1)**kronecker_delta(idx1,idx2)*(u_i*u_j-kronecker_delta(i,j))/l_u
    if m == 3:
        a,b,c = internal_coordinate
        u = get_vector(xyz_coordinates,a,b,False)
        l_u = np.linalg.norm(u)
        u /= l_u # Normalized u
        v = get_vector(xyz_coordinates,c,b,False)
        l_v = np.linalg.norm(v)
        v /= l_v # Normalized v
        cos_q = np.sum(u*v)
        sin_q = np.sqrt(1-cos_q**2)
        vector_dict = get_wilsonB_vector(xyz_coordinates,internal_coordinate)
        # m: a, o: b, n: c
        for indices in scanning_indices:
            x,y = indices
            idx1 = internal_coordinate[x]
            idx2 = internal_coordinate[y]
            # a: idx1, b: idx2
            for i in range(3):
                for j in range(3):
                    u_i = u[i]
                    u_j = u[j]
                    v_i = v[j]
                    v_j = u[j]
                    first_term = (u_i*v_j+u_j*v_i-3*u_i*u_j*cos_q+kronecker_delta(i,j)*cos_q)/(l_u**2*sin_q)
                    first_term *= sign_factor(idx1,a,b) * sign_factor(idx2,a,b)
                    second_term = (v_i*u_j+v_j*u_i-3*v_i*v_j*cos_q+kronecker_delta(i,j)*cos_q)/(l_v**2*sin_q)
                    second_term *= sign_factor(idx1,c,b) * sign_factor(idx2,c,b)
                    third_term = (u_i*u_j+v_j*v_i-u_i*v_j*cos_q-kronecker_delta(i,j))/(l_u*l_v*sin_q)
                    third_term *= sign_factor(idx1,a,b) * sign_factor(idx2,c,b)
                    fourth_term = (v_i*v_j+u_j*u_i-v_i*u_j*cos_q-kronecker_delta(i,j))/(l_u*l_v*sin_q)
                    fourth_term *= sign_factor(idx1,c,b) * sign_factor(idx2,a,b)
                    fifth_term = cos_q/sin_q * vector_dict[idx1][i] * vector_dict[idx2][j]
                    tensor[x][y][i][j] = first_term + second_term + third_term + fourth_term - fifth_term

    elif m == 4:
        a,b,c,d = internal_coordinate
        u = get_vector(xyz_coordinates,a,b,False)
        l_u = np.linalg.norm(u)
        u /= l_u # Normalized u
        v = get_vector(xyz_coordinates,d,c,False)
        l_v = np.linalg.norm(v)
        v /= l_v # Normalized v
        w = get_vector(xyz_coordinates,c,b,False)
        l_w = np.linalg.norm(w)
        w /= l_w # Normalized w
        uw = get_cross_vector(u,w,False)
        vw = get_cross_vector(v,w,False)
        cos_phi_u = np.sum(u*w)
        sin_phi_u = np.sqrt(1-cos_phi_u**2)
        cos_phi_v = np.sum(v*w)
        sin_phi_v = np.sqrt(1-cos_phi_v**2)
        cos_q = np.sum(uw*vw)/(sin_phi_u*sin_phi_v)
        # m: a, o: b, p: c, n: d
        for indices in scanning_indices:
            x,y = indices
            idx1 = internal_coordinate[x]
            idx2 = internal_coordinate[y]
            # a: idx1, b: idx2
            for i in range(3):
                for j in range(3):
                    u_i = u[i]
                    u_j = u[j]
                    v_i = v[j]
                    v_j = u[j]
                    w_i = w[i]
                    w_j = w[j]
                    first_term = (uw[i]*(cos_phi_u*w_j - u_j))/(l_u**2*sin_phi_u**4)
                    first_term *= sign_factor(idx1,a,b)*sign_factor(idx2,a,b)
                    second_term = (uw[j]*(cos_phi_u*w_i - u_i))/(l_u**2*sin_phi_u**4)
                    second_term *= sign_factor(idx1,a,b)*sign_factor(idx2,a,b)
                    third_term = (vw[i]*(w_j*cos_phi_v - v_j))/(l_v**2*sin_phi_v**4)
                    third_term *= sign_factor(idx1,d,c) * sign_factor(idx2,d,c)
                    fourth_term = (vw[j]*(w_i*cos_phi_v - v_i))/(l_v**2*sin_phi_v**4)
                    fourth_term *= sign_factor(idx1,d,c) * sign_factor(idx2,d,c)
                    fifth_term = (uw[i]*(w_j-2*u_j*cos_phi_u+w_j*cos_phi_u**2)) /(2*l_u*l_w*sin_phi_u**4)
                    fifth_term *= (sign_factor(idx1,a,b)*sign_factor(idx2,b,c) + sign_factor(idx1,c,b)*sign_factor(idx2,b,a))
                    sixth_term = (uw[j]*(w_i-2*u_i*cos_phi_u+w_i*cos_phi_u**2)) /(2*l_u*l_w*sin_phi_u**4)
                    sixth_term *= (sign_factor(idx1,a,b)*sign_factor(idx2,b,c) + sign_factor(idx1,c,b)*sign_factor(idx2,b,a))
                    seventh_term = (vw[i]*(w_j+2*u_j*cos_phi_v+w_j*cos_phi_v**2))/(2*l_v*l_w*sin_phi_v**4) 
                    seventh_term *= (sign_factor(idx1,d,c)*sign_factor(idx2,c,b)+sign_factor(idx1,c,b)*sign_factor(idx2,d,c))
                    eighth_term = (vw[j]*(w_i+2*u_i*cos_phi_v+w_i*cos_phi_v**2))/(2*l_v*l_w*sin_phi_v**4) 
                    eighth_term *= (sign_factor(idx1,d,c)*sign_factor(idx2,c,b)+sign_factor(idx1,c,b)*sign_factor(idx2,d,c))
                    ninth_term = (uw[i]*(u_j+u_j*cos_phi_u**2-3*w_j*cos_phi_u+w_j*cos_phi_u**3))/(2*l_w**2*sin_phi_u**4)
                    ninth_term *= sign_factor(idx1,b,c)*sign_factor(idx2,c,b)
                    tenth_term = (uw[j]*(u_i+u_i*cos_phi_u**2-3*w_i*cos_phi_u+w_i*cos_phi_u**3))/(2*l_w**2*sin_phi_u**4)
                    tenth_term *= sign_factor(idx1,b,c)*sign_factor(idx2,c,b)
                    eleventh_term = (vw[i]*(v_j+v_j*cos_phi_v**2+3*w_j*cos_phi_v-w_j*cos_phi_v**3))/(2*l_w**2*sin_phi_v**4)
                    eleventh_term *= sign_factor(idx1,b,c)*sign_factor(idx2,b,c)
                    twelvth_term = (vw[j]*(v_i+v_i*cos_phi_v**2+3*w_i*cos_phi_v-w_i*cos_phi_v**3))/(2*l_w**2*sin_phi_v**4)
                    twelvth_term *= sign_factor(idx1,b,c)*sign_factor(idx2,b,c)
                    if i != j:
                        k = list(set([0,1,2])-set([i,j]))[0]
                        w_k = w[k]
                        u_k = u[k]
                        v_k = v[k]
                        thirteenth_term = (j-i)*(-1/2)**(abs(i-j)) * (w_k*cos_phi_u-u_k)/(l_u*l_w*sin_phi_u)
                        thirteenth_term *= (1-kronecker_delta(idx1,idx2))*(sign_factor(idx1,a,b)*sign_factor(idx2,b,c)+sign_factor(idx1,c,b)*sign_factor(idx2,b,a))
                        fourteenth_term = (j-i)*(-1/2)**(abs(i-j))*(w_k*cos_phi_v-v_k)/(l_v*-l_w*sin_phi_v)
                        fourteenth_term *= (1-kronecker_delta(idx1,idx2))*(sign_factor(idx1,d,b)*sign_factor(idx2,b,c)+sign_factor(idx1,c,b)*sign_factor(idx2,b,a))
                    else:
                        thirteenth_term = 0
                        fourteenth_term = 0
                    tensor[x][y][i][j] = first_term + second_term + third_term + fourth_term + fifth_term + sixth_term + seventh_term + eighth_term + ninth_term + tenth_term + eleventh_term + twelvth_term + thirteenth_term + fourteenth_term

    return tensor, scanning_indices


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

def get_wilsonB_derivative_matrix(xyz_coordinates,internal_coordinates):
    n = len(xyz_coordinates)
    m = len(internal_coordinates)
    B_derivative_matrix = np.zeros((m,3*n,3*n))
    for row_idx, internal_coordinate in enumerate(internal_coordinates):
        tensor,scanning_indices = get_wilsonB_derivative_sub_matrix(xyz_coordinates,internal_coordinate) # tensor: m x m x 3 x 3
        for indices in scanning_indices:
            x,y = indices
            a = internal_coordinate[x]
            b = internal_coordinate[y]
            B_derivative_matrix[row_idx,3*a:3*a+3,3*b:3*b+3] = tensor[x,y,:,:]
    return B_derivative_matrix

def get_gradient_q(g_x,xyz_coordinates,internal_coordinates):
    B = get_wilsonB_matrix(xyz_coordinates,internal_coordinates) # M x 3N, B.T: 3N x M
    BT_inverse = np.linalg.pinv(B.T)
    return BT_inverse@g_x # M x 1

def get_hessian_q(g_x,h_x,xyz_coordinates,internal_coordinates):
    #print ('bf',internal_coordinates)
    B = get_wilsonB_matrix(xyz_coordinates,internal_coordinates) # M x 3N
    #print ('af',internal_coordinates)
    B_prime = get_wilsonB_derivative_matrix(xyz_coordinates,internal_coordinates) # M x 3N x 3N
    BT_inverse = np.linalg.pinv(B.T)
    B_inverse = np.linalg.pinv(B)
    g_q = BT_inverse@g_x # M x 1
    #print (g_q.shape,B_prime.shape)
    K = np.squeeze(np.einsum('il,ijk->ljk',g_q,B_prime),axis=0) # M x 1, M x 3N x 3N -> 1 x 3N x 3N (einsum) -> 3N x 3N (squeeze)
    #print (K.shape)
    #print (h_x.shape)
    #print ('K',np.sum(np.abs(K)))
    h_q = BT_inverse@(h_x - K)@B_inverse
    #print (h_q.shape)
    return h_q


def solve_equation(B,delta_q):
    B_inverse = np.linalg.pinv(B) # N x M
    delta_x = B_inverse@delta_q
    return delta_x


def get_xyz_updates(xyz_coordinates,internal_coordinates,delta_q): # Note that here, delta_qs[internal_coordinate] = delta_q
    B = get_wilsonB_matrix(xyz_coordinates,internal_coordinates)
    delta_x = solve_equation(B,delta_q) # 3N x 1 -> N x 3
    n = len(xyz_coordinates)
    delta_x = np.reshape(delta_x,(n,3))
    return delta_x


def update_geometry(xyz_coordinates,q_updates,criteria = 1e-6,max_iteration = 25):
    internal_coordinates = list(q_updates.keys())
    delta_q = np.array(list(q_updates.values())).T # M x 1
    q = get_q(xyz_coordinates,internal_coordinates) # q_n or q_0
    #print ('q: ',q)
    #print ('internals: ',internal_coordinates)
    target_q = q + delta_q # q
    original_coordinates = np.copy(xyz_coordinates)
    for i in range(max_iteration):
        delta_x = get_xyz_updates(xyz_coordinates,internal_coordinates,delta_q)
        xyz_coordinates += delta_x
        q = get_q(xyz_coordinates,internal_coordinates)
        delta_q = target_q - q # dq = q - q_n
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
    from pyMCD.utils import process
    molecule = chem.Molecule('R.xyz')
    coordinate_list = molecule.get_coordinate_list()
    #print (internal_coordinates)
    #exit()
    delta_qs = dict()
    #delta_qs[(0,1)] = 0.0
    #delta_qs[(1,2)] = 0.0
    #delta_qs[(2,3)] = 0.0
    internal_coordinate = (0,1,2)
    delta_qs[internal_coordinate] = np.pi/4
    internal_coordinate = (1,2,3)
    delta_qs[internal_coordinate] = np.pi/4
    internal_coordinate = (0,1)
    delta_qs[internal_coordinate] = 0.0
    internal_coordinate = (2,3)
    delta_qs[internal_coordinate] = 0.0
    internal_coordinate = (1,2)
    delta_qs[internal_coordinate] = 0.0
    internal_coordinate = (0,1,2,3)
    delta_qs[internal_coordinate] = -np.pi*3/4

    #delta_qs[(0,1,2,3)] = np.pi/3
    #delta_qs[(1,2,3)] = 0.0
    #delta_qs[(0,1,2,3)] = np.pi*2/3
    internal_coordinates = list(delta_qs.keys())
    #print (get_q(coordinate_list,internal_coordinates))
    matrix = get_wilsonB_derivative_matrix(coordinate_list,internal_coordinates)
    #print (matrix)
    #print (matrix.shape)
    
    update_geometry(coordinate_list,delta_qs)
    process.locate_molecule(molecule,coordinate_list)
    molecule.print_coordinate_list()

    #for i,coordinate in enumerate(coordinate_list):
    #    print (molecule.atom_list[i].get_element(), coordinate[0],coordinate[1],coordinate[2])
    #print (get_internal_coordinate_info(coordinate_list,internal_coordinates))
