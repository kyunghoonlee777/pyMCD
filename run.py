
### python provided modules ###
import time
import os
import sys
from copy import deepcopy
import pickle
import datetime
import argparse

### extra common libraries ###
import numpy as np

### ace-reaction libraries ###
from pyMCD import chem
from pyMCD import mcd
from pyMCD.utils import process


def read_reactant(directory):
    f = open(os.path.join(directory,'R.com'))
    state_info = f.readline().strip().split(' ')
    chg,multiplicity = int(state_info[0]), int(state_info[1])
    atom_list = []
    atom_info = f.readline()
    while atom_info.strip() != '':
        atom_info = atom_info.strip().split()
        print (atom_info)
        atom_type = atom_info[0]
        x = float(atom_info[1])
        y = float(atom_info[2])
        z = float(atom_info[3])
        atom = chem.Atom(atom_type)
        atom.x = x
        atom.y = y
        atom.z = z
        atom_list.append(atom)
        try:
            atom_info = f.readline()
            if atom_info.strip() == '':
                break
        except:
            break
    f.close()
    reactant = chem.Molecule()
    reactant.atom_list = atom_list
    reactant.adj_matrix = process.get_adj_matrix_from_distance(reactant)
    return reactant,chg,multiplicity

def read_bond_info(directory):
    constraints = dict()
    num_steps = dict()
    formed_bonds = []
    broken_bonds = []
    with open(os.path.join(directory,'bond_info')) as f:
        for line in f:
            info = line.strip().split() #0: start, 1: end, 2: target length, 3: Num steps
            constraint = tuple([int(idx) - 1 for idx in info[:-2]])
            target_value = float(info[-2])
            if len(constraint) > 2:
                target_value *= np.pi/180
            constraints[constraint] = target_value
            num_steps[constraint] = int(info[-1])
    print (constraints)
    return constraints, num_steps


def change_option(args):
    option_directory = os.path.join(args.input_directory,'option')
    if os.path.exists(option_directory):
        with open(option_directory) as f:
            for line in f:
                words = line.strip().split('=')
                attribute = words[0]
                value = words[1]
                if attribute == 'working_directory':
                    args.working_directory = value
                if attribute == 'num_relaxation':
                    args.num_relaxation = int(value)
                if attribute == 'scale':
                    args.scale = float(value)
                if attribute == 'unit':
                    args.unit = value
    else:
        print ('option directory is not found! Default parameters are used!!!')

def get_calculator(args):
    calculator_name = args.calculator.lower()
    if calculator_name == 'gaussian':
        from pyMCD.Calculator import gaussian
        calculator = gaussian.GaussianMCD('g16')
    elif calculator_name == 'orca':
        from pyMCD.Calculator import orca
        calculator = orca.Orca()
    calculator.load_content(os.path.join(args.input_directory,'qc_input'))
    basis_file = os.path.join(args.input_directory,'basis') # For Effective Core Potential
    if os.path.exists(basis_file):
        calculator.load_basis(basis_file)
    return calculator


def generate_path():
    import datetime
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_directory','-id',type=str,help='directory of inputs')
    parser.add_argument('--save_directory','-sd',type=str,help='directory for saving outputs',default=None)
    parser.add_argument('--working_directory','-wd',type=str,help='scratch directory of QC programs',default='')
    parser.add_argument('-num_relaxation',type=int,help='Num relaxation for input',default=5)
    parser.add_argument('-scale',type=float,help='Maxmial displacement',default=1.0)
    parser.add_argument('--calculator','-c',type=str,help='Name of Quantum Calculation software',default='gaussian')
    parser.add_argument('--unit','-u',type=str,help='unit',default='Hartree')

    args = parser.parse_args()

    input_directory = args.input_directory
    save_directory = args.save_directory
    if save_directory is None:
        save_directory = input_directory
    reactant, chg, multiplicity = read_reactant(input_directory) # Read geometry of reactant
    constraints, num_steps = read_bond_info(input_directory) # bond info
    change_option(args) # Read option file and change values in args
    calculator = get_calculator(args) # Make calculator, you can use your own calculator!
    scanner = mcd.MCD(num_relaxation = args.num_relaxation,calculator=calculator)
    scanner.scale = args.scale
    if args.working_directory != '':
        scanner.change_working_directory(args.working_directory)
    scanner.change_energy_unit(args.unit)
    print (chg,multiplicity,len(reactant.atom_list))
    scanner.log_directory = save_directory
    # Also, write constraints information (TODO)
    scanner.consider_redundant = False
    pathway = scanner.scan(reactant,constraints,num_steps,chg = chg, multiplicity = multiplicity)

if __name__ == '__main__':
    generate_path()
    print ('hi')
