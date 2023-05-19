'''
--chem.py--
Declaring classes.
(Intermediates, molecules, atoms, ...)
Some parts can be replaced by ASE formats, but now we are using these settings for convenient customization.
'''

import numpy as np

from pyMCD.utils import ic
from pyMCD.utils import process


class Atom:
    """
    :class Atom:
        class Atom mainly contains atomic_number, element, and x,
        Other attributes are not widely used
        molecule_index shows on which molecule that this atom is contained in. For example, if we consider Intermediate C.C, 
        molecule_index can have either 0 or 1, and every atom within atom_list of Intermediate can be assigned by checking
        which molecule does the atom belong to.

    :param data(str or integer):
        Data is either symbol that represents element or atomic number
    
    """
    global periodic_table
    periodic_table = ['H','He','Li','Be','B','C','N','O','F','Ne',\
    'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn',\
    'Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr',\
    'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba',\
    'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',\
    'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn']

    def __init__(self,data = None):
        self.atomic_number = None
        self.element = None
        self.x = 0.00
        self.y = 0.00
        self.z = 0.00
        molecule_index = 0
        self.is_divalent_hydrogen = False
        if data is not None:
            if type(data) == str:
                self.element = data
            else:
                self.atomic_number = data


    def set_atomic_number(self,atomic_number):
        self.atomic_number = atomic_number
        self.element = periodic_table[atomic_number-1]

    def set_element(self, element):
        """ Type of an atom. e.g. 'C', 'H', 'O', and so on.""" 
        self.element = element
        self.atomic_number = periodic_table.index(element) + 1

    def set_coordinate(self, position):
        """ Set Cartesian coordinates of an atom """ 
        dim = len(position)
        if dim == 2:
            x = position[0]
            y = position[1]
            z = 0
        elif dim == 3:
            x = position[0]
            y = position[1]
            z = position[2]
        self.x = x
        self.y = y
        self.z = z

    def get_atomic_number(self):
        """
        Returns the atomic number (number of protons) of a given atom.

        :param:

        :return integer:
            Directly the atomic number of a given atom is returned
        """
        if self.atomic_number == None:
            element = self.element
            if element == None:
                print ('atom is not specified!')
            if len(element)>1:
                end_part = element[1:]
                end_part = str.lower(end_part)
                element = element[0]+end_part
                self.element = element
            if element in periodic_table:
                index = periodic_table.index(element)
                self.atomic_number = index + 1
            else:
                print ('element',element)
                print ('modify periodic table!!!')
        return self.atomic_number

    def get_element(self):
        """
        Returns symbol of a given atom.

        :param:

        :return str:
            Directly the symbol of a given atom is returned
        """
        if self.element == None:
            atomic_number = self.atomic_number
            if atomic_number == None:
                print ('atom is not specified!')
            z = int(self.atomic_number)-1
            self.element = periodic_table[z]
        return self.element

    def get_coordinate(self):
        return np.array([self.x,self.y,self.z])           

    def get_mass(self):
        """
        Returns the exact value (float) of an atomic mass of a given atom.

        :param:

        :return float:
            Directly the atomic mass of a given atom is returned
        """
        element = self.get_element()
        a = str.lower(element)
        if a=='h': return 1.008
        elif a=='he': return 4.003
        elif a=='li': return 6.941
        elif a=='be': return 9.012
        elif a=='b': return 10.81
        elif a=='c': return 12.01
        elif a=='n': return 14.01
        elif a=='o': return 16.00
        elif a=='f': return 19.00
        elif a=='ne': return 20.18
        elif a=='na': return 22.99
        elif a=='mg': return 24.31
        elif a=='al': return 26.98
        elif a=='si': return 28.09
        elif a=='p': return 30.97
        elif a=='s': return 32.07
        elif a=='cl': return 35.45
        elif a=='ar': return 39.95
        elif a=='k': return 39.10
        elif a=='ca': return 40.08
        elif a=='au': return 196.97
        elif a=='co': return 58.9332
        elif a=='ni': return 58.6934
        elif a=='ti': return 47.8671
        elif a=='fe': return 55.845
        elif a=='br': return 79.904
        elif a=='rh': return 102.90550
        elif a=='pd': return 106.42
        elif a=='hf': return 178.49
        elif a=='i': return 126.90447
        else: return 0

    def get_radius(self):
        """
        Returns a radius information of a given atom. Reference is given here: Dalton Trans., 2008, 2832-2838
        
        :param:

        :return float:
            It directly returns the reference values
        """
        element = self.get_element()       
        a = str.lower(element)
        #reference : Dalton Trans., 2008, 2832-2838
        if a=='h': 
            return 0.31
        elif a=='li': 
            return 1.28
        elif a=='be': 
            return 0.96
        elif a=='b': 
            return 0.84
        elif a=='c': 
            return 0.76
        elif a=='n': 
            return 0.71
        elif a=='o': 
            return 0.66
        elif a=='f': 
            return 0.57
        elif a=='na': 
            return 1.66
        elif a=='mg': 
            return 1.41
        elif a=='al': 
            return 1.21
        elif a=='si': 
            return 1.11
        elif a=='p': 
            return 1.07
        elif a=='s': 
            return 1.05
        elif a=='cl': 
            return 1.02
        elif a=='ar': 
            return 0.76
        elif a=='k': 
            return 2.03
        elif a=='ca': 
            return 1.76
        elif a=='co': 
            return 1.38 #1/2*(lowspin+highspin)
        #elif a=='co': return 1.26 #lowspin
        #elif a=='co': return 1.50 #highspin
        elif a=='fe': 
            return 1.42 #1/2*(lowspin+highspin)
        elif a=='ni': 
            return 1.24
        #elif a=='cr': return 1.39
        elif a=='ti': 
            return 1.60
        elif a=='br': 
            return 1.20
        elif a=='rh': 
            return 1.42
        elif a=='pd': 
            return 1.39
        elif a=='i': 
            return 1.39
        elif a=='hf': 
            return 1.75
        else: 
            return 0

        #reference : J. Chem. Phys. 41, 3199 (1964)
        '''
        if a=='h': return 0.25
        elif a=='li': return 1.45
        elif a=='be': return 1.05
        elif a=='b': return 0.85
        elif a=='c': return 0.70
        elif a=='n': return 0.65
        elif a=='o': return 0.60
        elif a=='f': return 0.50
        elif a=='na': return 1.80
        elif a=='mg': return 1.50
        elif a=='al': return 1.25
        elif a=='si': return 1.10
        elif a=='p': return 1.00
        elif a=='s': return 1.00
        elif a=='cl': return 1.00
        elif a=='ar': return 0.71
        elif a=='k': return 2.20
        elif a=='ca': return 1.80
        elif a=='co': return 1.35
        else: return 0
        '''

    def get_period_group(self):
        """
        Returns a period,group information from a given atom. It finds period, group by identifying 
        electron configuration (orbital configuration)
        
        :param:

        :return period,group(int,int):
            Note that these values are actual period and group. If C atomm is given, it returns 2,4

        """
        atomic_number = self.get_atomic_number()
        num_of_electrons = atomic_number
        # Orbital: [n,l,num_of_electrons]
        sum_of_n_and_l=1
        orbital_configuration=[]
        while num_of_electrons>0:
            # Generate orbitals within sum_of_n_and_l=k
            # New orbitals are introduced for (k+1)/2
            maximal_l=int((sum_of_n_and_l-1)/2)
            for l in range(maximal_l,-1,-1):
                # Start with lowest l
                if num_of_electrons>4*l+2:
                    num_of_electrons-=(4*l+2)
                    orbital_configuration.append([sum_of_n_and_l-l,l,4*l+2])
                else:
                    orbital_configuration.append([sum_of_n_and_l-l,l,num_of_electrons])
                    num_of_electrons=0
                    break
            sum_of_n_and_l+=1
        # Get maximal n and l
        period=0
        for orbital in orbital_configuration:
            if orbital[0]>period:
                period = orbital[0]
        # If transition metal, we add 9, Sc has group 9 for else, we do not consider ...
        last_orbital = orbital_configuration[-1]
        if last_orbital[1]<2:
            group = 2 * last_orbital[1] ** 2 + last_orbital[2]
        else:
            group = 8 + last_orbital[2]
        return period,group

    def copy(self):
        new_atom = Atom()
        # Copy all attributes
        new_atom.atomic_number = self.atomic_number
        new_atom.element = self.element
        new_atom.x = self.x
        new_atom.y = self.y
        new_atom.z = self.z
        return new_atom 

    def is_same_atom(self,atom):
        """
        Returns whether the two atoms have same type by comparing atomic number

        :param atom(pyclass 'Atom'):
            Our defined class 'Atom'

        :return:
            True: Two atoms are same type
            False: Two atoms are different type
        """
        atomic_number = self.get_atomic_number()
        atomic_number_prime = atom.get_atomic_number()
        return atomic_number == atomic_number_prime

    def get_content(self,option='element',criteria = 1e-4):
        x = self.x
        y = self.y
        z = self.z
        if abs(x) < criteria:
            x = 0.00
        if abs(y) < criteria:
            y = 0.00
        if abs(z) < criteria:
            z = 0.00
        content = ' ' + str(x) + ' ' + str(y) + ' ' + str(z) + '\n'
        if option=='element':
            content = self.get_element() + content
        else:
            content = str(self.get_atomic_number()) + content
        return content

    def __eq__(self,atom):
        return self.is_same_atom(atom)
         
class Molecule:
    """
    :class Molecule:
        class Molecule mainly contains atom_list, atom_feature, chg, bo_matrix, adj_matrix, energy, smiles, c_eig_list
        atom_list is a list of atom (pyclass 'Atom')
        atom_feature is a dict, where features of atom such as formal charge, number of pi bonds, etc, are stored.
        Those information can be freely added by giving new dict
        c_eig_list can be used to identify whether the given two molecules are the same. This c_eig_list is invariant
        to the permutation of atom indexing, identification between any two generated molecules can be easily checked.

    :param data(str or xyz file or None):
        data should be feeded as either smiles(str), xyz file(file) or None
        * If smiles is used, rd_mol is generated using the given smiles and converted into ace_mol (pyclass 'Molecule)
        * If xyz file is used, directly the 3d geometry of the given molecule is generated. If you want to generate adj_matrix,
        bo_matrix, chg_list, etc, refer following method contained in pyclass 'Molecule' (charge of molecule should be given!!!)
        i. Generate adj_matrix by using 'get_adj_matrix_from_distance' stored in class 'Molecule'
        ii. Generate bo_matrix by using 'get_adj_matrix_from_adj_matrix' stored in class 'Molecule'
        iii. Then, using those bo_matrix, get chg_list by using 'get_chg_list_from_bo_matrix' stored n class 'Molecule'
        * If None is used, only blank virtual molecule is generated

    """
    def __init__(self,data = None):
        self.atom_list = []
        self.atom_feature = dict()
        self.adj_matrix = None
        self.bo_matrix = None
        self.chg = None
        self.energy = None        
        self.smiles = None
        self.multiplicity = None

        if data == None:
            pass
          
        elif (type(data) == str) and (data[-4:] == '.xyz'):
            # data is already opened file
            f = open(data,'r')
            atom_list=[]
            try:
                atom_num = int(f.readline())
            except:
                print ('Wrong format! Should start with number of atoms!')
            try:
                energy = float(f.readline())
                self.energy = energy
            except:
                self.energy = None
            for i in range(atom_num):
                try:
                    content = f.readline().strip()
                    atom_line = content.split()
                    #atomic_number = int(atom_line[0])
                    element_symbol = atom_line[0]
                    x = float(atom_line[1]) 
                    y = float(atom_line[2]) 
                    z = float(atom_line[3])
                    new_atom = Atom(element_symbol)
                    new_atom.x = x
                    new_atom.y = y
                    new_atom.z = z
                    atom_list.append(new_atom)
                except:
                    print ('Error found in:',content)
                    print ('Check the file again:',data)
            self.atom_list = atom_list
            # At least make adjacency
            #self.adj_matrix = process.get_adj_matrix_from_distance(self)
              
    def get_chg(self):
        if self.chg is None:
            chg_list = self.get_chg_list()
            if chg_list is None:
                return None
            else:
                return np.sum(chg_list)
        else:
            return self.chg

    def get_multiplicity(self):
        try:
            e_list = molecule.get_num_of_lone_pair_list()
            num_of_unpaired_e = len(np.where((2*e_list) % 2 == 1)[0])    
            multiplicity = num_of_unpaired_e + 1
            return multiplicity
        except:
            return None

    def get_adj_matrix(self):
        if self.adj_matrix is not None:
            return self.adj_matrix
        if self.bo_matrix is not None:
            adj_matrix = np.where(self.bo_matrix>0,1,0)
            return adj_matrix
        return None
    
    def get_bo_matrix(self):            
        return self.bo_matrix

    def copy(self):
        new_molecule = Molecule()
        atom_list = self.atom_list
        # First copy atoms
        new_atom_list = []
        for atom in atom_list:
            new_atom_list.append(atom.copy())
        new_molecule.atom_list = new_atom_list
        # Copy connectivity information
        bo_matrix = self.get_bo_matrix()
        if bo_matrix is not None:
            new_molecule.bo_matrix = np.copy(bo_matrix)
        else:
            adj_matrix = self.get_adj_matrix()
            if adj_matrix is not None:
                new_molecule.adj_matrix = np.copy(adj_matrix)
            else:
                print ('Warning: Connectivity information is not included in the molecule!!!')
        return new_molecule 

    def get_z_list(self):
        """
        Returns atomic number list of a given molecule
         
        :param:

        :return z_list(pyclass 'numpy.ndarray'):
            list of atomic number
        """
        atom_feature = self.atom_feature
        if atom_feature is not None and 'atomic number' in atom_feature:
            return atom_feature['atomic number']
        else:
            z_list = list(map(lambda x:x.get_atomic_number(),self.atom_list))
            return np.array(z_list)

    def get_element_list(self):
        """
        Returns element list of a given molecule
         
        :param:

        :return element_list(list of string):
            list of element written as capital letters
        """
        atom_feature = self.atom_feature
        if atom_feature is not None and 'element' in atom_feature:
            return atom_feature['element']
        else:
            element_list = list(map(lambda x:x.get_element(),self.atom_list))
            return element_list


    def get_chg_list(self):
        """
        Returns chg list of a given molecule
         
        :param:

        :return chg_list(pyclass 'numpy.ndarray'):
            list of formal charge
        """
        atom_feature = self.atom_feature
        try:
            return atom_feature['chg']
        except:
            print ('charge lists are not prepared!!!')
            return None

    def get_radius_list(self):
        """
        Returns radius list of a given molecule
        Ex. For CH4 with atom order [C,H,H,H,H], then radius_list is given as
        [0.8,0.5,0.5,0.5,0.5], if radius of C and H are given as 0.8 and 0.5 (unit is Angstrom)

        :param:

        :return radius_list(list of float):
            list of radius of each atom
        """
        atom_list = self.atom_list
        atom_feature = self.atom_feature
        if 'radius' in atom_feature:
            return atom_feature['radius']
        radius_list = []
        for atom in atom_list:
            radius_list.append(atom.get_radius())
        return radius_list

    def get_bond_list(self,contain_bond_order = True):
        """
        Returns the total bond list as list of tuples
        For example, if CH4 is given with atom order [C,H,H,H,H], if contain_bond_order = False, the output is given as
        [(0,1),(0,2),(0,3),(0,4)]
        if contain_bond_order = True, the output is given as
        [(0,1,1),(0,2,1),(0,3,1),(0,4,1)]

        :param contain_bond_order(boolean):
            If contain_bond_order is False, it only returns bonds represented as (i,j) between atoms within the given intermediate.
            If contain_bond_order is True, it returns bonds with bond order included (i,j,bond_order), where bond_order can only have 1,2,3.
            Therefore, a given molecule(self) should be kekulized.

        :return bond_list(either list of tuple with size 2 or 3):
            bond_list
        """
        atom_list = self.atom_list
        n = len(atom_list)
        check_matrix = self.bo_matrix
        total_bond_list = []
        if contain_bond_order:
            check_matrix = self.get_bo_matrix()
            if check_matrix is None:
                print ('we cannot give bond order!!!')
                print ('We will automatically give only bond list!')
                contain_bond_order = False
                check_matrix = self.get_adj_matrix()
                if check_matrix is None:
                    print ('matrix',self.atom_list)
                    print ('hahahahahaha',check_matrix)
                    print ('Give connectivity! We cannot find the bond!')
                    return None
        if contain_bond_order:            
            bond_type = [1,2,3]
            check_matrix = self.get_bo_matrix()
        else:
            bond_type = [1]
            check_matrix = self.get_adj_matrix()
        #check_matrix = self.adj_matrix
        for bond_order in bond_type:
            bond_list = np.where(check_matrix == bond_order)
            bond_list = np.stack(bond_list,axis = 1)
            for array in bond_list:
                if array[0] < array[1]:
                    if contain_bond_order:
                        bond_tuple = (int(array[0]),int(array[1]),int(bond_order))
                    else:
                        bond_tuple = (int(array[0]),int(array[1]))
                    total_bond_list.append(bond_tuple) 
        return total_bond_list 

    def get_minimal_data(self):
        data = dict()
        data['z'] = self.get_z_list()
        data['adj'] = self.get_adj_matrix()
        data['bo'] = self.get_bo_matrix()
        data['chg'] = self.get_chg()
        data['atom chg'] = self.get_chg_list()
        data['coords'] = self.get_coordinate_list()
        return data

    def get_coordinate_list(self):
        """
        Returns 3d coordinate of a given molecule
        
        :param:

        :return coordinate_list(list(size n) of tuple of float(size 3)):
            
        """
        coordinate_list = []
        atom_list = self.atom_list
        if 'coords' in self.atom_feature:
            return self.atom_feature['coords']
        for atom in atom_list:
            coordinate_list.append([atom.x,atom.y,atom.z])
        return np.array(coordinate_list)


    def print_coordinate_list(self,option='element'):
        coordinate_list = self.get_coordinate_list()
        atom_list = self.atom_list
        n = len(atom_list)
        for i in range(n):
            coordinate = coordinate_list[i]
            element = atom_list[i].get_element()
            if option == 'number':
                element = atom_list[i].get_atomic_number()
            print_x = coordinate[0]
            print_y = coordinate[1]
            print_z = coordinate[2]
            if abs(print_x) < 0.0001:
                print_x = 0.00
            if abs(print_y) < 0.0001:
                print_y = 0.00
            if abs(print_z) < 0.0001:
                print_z = 0.00
            print (element, print_x,print_y,print_z)

    def write_geometry(self,file_directory, option='element',criteria = 1e-4):
        """
        Writes xyz file that contains the 3d molecular geometry
    
        :param file_directory(str):
            Directory for saving 3d geometry xyz file

        :return:

        """
        atom_list = self.atom_list
        n = len(atom_list)
        f = open(file_directory, 'w')
        if True: # If inappropriate geometry condition is determined, it will be added
            content = str(n)+'\n'            
            if self.energy is not None:
                content = content + str(self.energy)+'\n'
            else:
                content = content + '\n'
            f.write(content)
            for atom in atom_list:
                f.write(atom.get_content(option,criteria))
            f.close()
        else:
            print ('Wrong geometry!!!')        

    def get_normal_vector(self,idx1,idx2,idx3):
        """
        For given three atom indices, it returns normal vector that is perpendicular to
        the plane generated by three indices

        :param idx1,idx2,idx3(int):
            Indices for selected three atoms

        :return normal_vector(pyclass 'numpy.ndarray' with length 3)
            normal vector
        """
        vector1 = self.get_vector_between_atoms(idx3,idx1)
        vector2 = self.get_vector_between_atoms(idx3,idx2)
        cross_vector = np.cross(vector1,vector2)
        norm = np.linalg.norm(cross_vector)
        if norm == 0: 
            return cross_vector
        else:
            return cross_vector/norm
   
    def get_vector_between_atoms(self,idx1,idx2,normalize = True):
        """
        For given two atom indices, it returns difference vector between two atoms. This function can be used to 
        move or rotate a given molecule
        
        :param idx1,idx2(int):
            Indices for selected two atoms 

        :return difference_vector(pyclass 'numpy.ndarray' with length 3)

        """
        atom_list = self.atom_list
        atom_coord1 = atom_list[idx1].get_coordinate()
        atom_coord2 = atom_list[idx2].get_coordinate()
        vector = atom_coord2 - atom_coord1
        if normalize:
            norm = np.linalg.norm(vector)
            if norm < 0.0001:
                print ('zero vector is found ...')
                return vector
            else:
                #print ('finalvector',vector,vector/norm)
                return vector/norm
        return vector

    def get_internal_coordinate(self,indices,unit='degree'):
        if len(indices) == 2:
            idx1,idx2 = indices
            return self.get_distance_between_atoms(idx1,idx2)
        elif len(indices) == 3:
            idx1, idx2, idx3 = indices
            return self.get_angle_between_atoms(idx1,idx2,idx3,unit)
        elif len(indices) == 4:
            idx1, idx2, idx3, idx4 = indices
            return self.get_dihedral_angle_between_atoms(idx1,idx2,idx3,idx4,unit)
        else:
            print (f'Wrong coordinate (={indices}) given!')
            return None


    def get_distance_between_atoms(self,idx1,idx2):
        """
        Returns the distance between chosen two atoms

        :param idx1,idx2(int):
            indices of chosen two atoms. 

        :return distance(float):
            Distance between selected two atoms

        """
        coordinate_list = self.get_coordinate_list()
        return ic.get_distance(coordinate_list,idx1,idx2)

    def get_angle_between_atoms(self,idx1,idx2,idx3,unit='rad'):
        """
        Returns the distance between chosen two atoms

        :param idx1,idx2(int):
            indices of chosen two atoms. 

        :return distance(float):
            Distance between selected two atoms

        """
        coordinate_list = self.get_coordinate_list()
        angle = ic.get_angle(coordinate_list,idx1,idx2,idx3)
        if unit == 'degree':
            angle *= 180/np.pi
        return angle


    def get_dihedral_angle_between_atoms(self,idx1,idx2,idx3,idx4,unit='rad'):
        coordinate_list = self.get_coordinate_list()
        angle = ic.get_dihedral_angle(coordinate_list,idx1,idx2,idx3,idx4)
        if unit == 'degree':
            angle *= 180/np.pi
        return angle

