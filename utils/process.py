
"""
---generate_molecule.py---
Generate 2D graph molecule in ACE-Reaction format from xyz files or SMILES
"""

import os 
import subprocess
import copy
import numpy as np 
import random

from rdkit import Chem
import rdkit.Chem.rdchem as rdchem

from pyMCD import chem

from scipy import spatial


def locate_molecule(ace_molecule,coordinate_list,update = False):
    """ Locates atoms according to coordinate_list, be cautious on ordering of atoms
    Args:
        |  ace_molecule (<class 'Molecule'>): class instance of Molecule
        |  coordinate_list (list of list (x,y,z)): list of 3d coordinate, where each x,y,z are float type
    Returns:
        |  No return, it direcly modifies the geometry of a given molecule
    """
    atom_list = ace_molecule.atom_list
    for i in range(len(atom_list)):
        atom = atom_list[i]
        atom.set_coordinate(coordinate_list[i])
    if update:
        ace_molecule.set_adj_matrix(None)

def translate_molecule(ace_molecule,vector):
    atom_list = ace_molecule.atom_list
    for atom in atom_list:
        translate_atom(atom,vector)


def locate_atom(atom,coordinate):
    """ Locates a single atom to input 'coordinate'
    Args:
        |  atom (<class 'Atom'>): class instance of Atom
        |  coordinate (list of float): coordinate in form of [x,y,z]
    Returns:
        |  No return, it directly locates the atom to the given coordinate
    """
    atom.x = coordinate[0]
    atom.y = coordinate[1]
    atom.z = coordinate[2]

def translate_atom(atom,vector):
    atom.x += vector[0]
    atom.y += vector[1]
    atom.z += vector[2]

def get_ace_mol_from_rd_mol(rd_molecule,add_hydrogen = True,include_stereo = False):
    """ It converts rd_molecule type info ace_molecule type
    Args:
        |  rd_molecule (<class 'rdkit.Molecule>')
    Returns:
        |  ace_molecule(<class Molecule>)
    """
    # Kekulize molecule
    from rdkit import Chem
    rd_molecule_copy = copy.deepcopy(rd_molecule)
    try:
        if add_hydrogen:
            rd_molecule = Chem.AddHs(rd_molecule) 
        Chem.rdmolops.Kekulize(rd_molecule)
    except:
        rd_molecule = rd_molecule_copy
    bond_types = {Chem.BondType.SINGLE:1, Chem.BondType.DOUBLE:2, Chem.BondType.TRIPLE:3} 
    n = rd_molecule.GetNumAtoms()
    atom_list = []
    chg_list = []
    atom_feature = dict()
    # Make atom_list
    for i in range(n):
        rd_atom = rd_molecule.GetAtomWithIdx(i)
        ace_atom = chem.Atom()
        chg_list.append(rd_atom.GetFormalCharge())
        '''
        position = rd_molecule.GetAtomPosition(i)
        if position!=None:
            ace_atom.x = position[0]
            ace_atom.y = position[1]
            ace_atom.z = position[2]
        '''
        ace_atom.atomic_number = rd_atom.GetAtomicNum()
        atom_list.append(ace_atom)
    atom_feature['chg'] = np.array(chg_list)
    # Make bond order matrix
    bonds = rd_molecule.GetBonds()
    bo_matrix = np.zeros((n,n))
    for bond in bonds:
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        bond_order = bond_types[bond.GetBondType()]
        bo_matrix[begin][end] = bo_matrix[end][begin] = bond_order
    ace_molecule = chem.Molecule()
    ace_molecule.atom_list = atom_list
    ace_molecule.bo_matrix = bo_matrix
    ace_molecule.atom_feature = atom_feature
    return ace_molecule

def get_adj_matrix_from_distance(molecule,coeff = 1.10):
    """
    Returns adj_matrix from 3d coordinate of given molecule
    It recognizes bond between two atoms, if the sum of radius * coeff is less than distance between two atoms

    :param coeff(float):
        criteria for recognizing bond. If criteria gets higher, more and more bonds are generated between atoms, 
        since criteria distance for bond distance gets higher.
        Appropriate criteria value is between 0.8 ~ 1.3, here we set default value as 1.10
    
    :return adj(pyclass 'numpy.ndarray'):
        connectivity matrix between atoms
         
    """
    atom_list = molecule.atom_list
    n = len(atom_list)
    radius_list = molecule.get_radius_list()
    radius_matrix_flatten = np.repeat(radius_list,n)
    radius_matrix = radius_matrix_flatten.reshape((n,n))
    radius_matrix_transpose = np.transpose(radius_matrix)
    criteria_matrix = coeff * (radius_matrix + radius_matrix_transpose)
    coordinate_list = molecule.get_coordinate_list()
    distance_matrix = spatial.distance_matrix(coordinate_list,coordinate_list)
    adj = np.where(distance_matrix<criteria_matrix,1,0)
    np.fill_diagonal(adj,0)
    return adj 


