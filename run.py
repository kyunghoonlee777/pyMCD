
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
from pyMCD import mcd_original
from pyMCD import mcd
from pyMCD.utils import process


def read_reactant(directory):
    try:
        f = open(os.path.join(directory,'R.com'))
    except:
        print ('Cannot find \'R.com\' file! Recheck your input !!!')
        exit()
    state_info = f.readline().strip().split(' ')
    chg,multiplicity = int(state_info[0]), int(state_info[1])
    atom_list = []
    atom_info = f.readline()
    while atom_info.strip() != '':
        atom_info = atom_info.strip().split()
        #print (atom_info)
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
    reactant.chg = chg
    reactant.multiplicity = multiplicity
    return reactant

def read_bond_info(directory):
    constraints = dict()
    num_steps = dict()
    formed_bonds = []
    broken_bonds = []
    try:
        f =  open(os.path.join(directory,'coordinates'))
    except:
        print ('Cannot find \'coordinates\' file! Recheck your input !!!')
        exit()
    for line in f:
        info = line.strip().split() #0: start, 1: end, 2: target length, 3: Num steps
        if len(info) > 0:
            try:
                constraint = tuple([int(idx) - 1 for idx in info[:-2]])
                target_value = float(info[-2])
            except:
                print ('Wrong coordinate found! Check coordinate file again !!!')
                print (f'Current line: {line.strip()}')
                print ('Input should be given as: [Coordinate] [target value] [num_step]')
                exit()
            if len(constraint)  < 2:
                print ('Wrong coordinate found! Check coordinate file again !!!')
                print (f'Current line: {line.strip()}')
                print ('Input should be given as: [Coordinate] [target value] [num_step]')
                exit()
            if len(constraint) > 2:
                target_value *= np.pi/180
            constraints[constraint] = target_value
            try:
                num_steps[constraint] = int(info[-1])
            except:
                print ('Num steps are not given! Check coordinate file again !!!')
                print ('Input should be given as: [Coordinate] [target value] [num_step]')
                exit()
    return constraints, num_steps


def change_option(args):
    correct = True
    try:
        option_directory = os.path.join(args.input_directory,'option')
    except:
        print ('Default option is used !!!')
        return
    if os.path.exists(option_directory):
        wrong_attributes = []
        with open(option_directory) as f:
            for line in f:
                words = line.strip().split('=')
                attribute = words[0]
                value = words[1]
                if attribute == 'working_directory':
                    args.working_directory = value
                elif attribute == 'num_relaxation':
                    try:
                        value = int(value)
                        if value <= 0:
                            print (f'Wrong num_relaxation (={value}) is given! Check the option file !!!')
                            print (f'The value must be positive integer !!!\n')
                            correct = False
                        else:
                            args.num_relaxation = value
                    except:
                        print (f'Wrong num_relaxation (={value}) is given! Check the option file !!!')
                        print (f'The value must be positive integer !!!\n')
                        correct = False
                elif attribute == 'step_size':
                    try:
                        value = float(value)
                        if value <= 0.0:
                            print (f'Wrong step_size (={value}) is given! Check the option file !!!')
                            print (f'The value must be positive !!!\n')
                            correct = False
                        else:
                            args.step_size = value
                    except:
                        print (f'Wrong step_size (={value}) is given! Check the option file !!!')
                        print (f'The value must be positive !!!\n')
                        correct = False
                elif attribute == 'unit':
                    if not value in ['eV','Hartree','kcal']:
                        print (f'Wrong unit (={value}) is given! Check the option file !!!')
                        print ('Only eV, Hartree, kcal are allowed options !!!\n')
                        correct = False
                    else:
                        args.unit = value
                elif attribute == 'calculator':
                   args.calculator = value
                elif attribute == 'command':
                    args.command = value
                elif attribute == 'use_hessian':
                    try:
                        value = int(value)
                        if not value in [0,1]:
                            print (f'Wrong use_hessian (={value}) is given! Check the option file !!!')
                            print ('Only 0 or 1 are possible. If the value is zero, hessian is not used. Otherwise, hessian is used ... Default value is 0 \n')
                            correct = False
                        else:
                            args.use_hessian = value
                    except:
                        print (f'Wrong use_hessian (={value}) is given! Check the option file !!!')
                        print ('Only 0 or 1 are possible. If the value is zero, hessian is not used. Otherwise, hessian is used ... Default value is 0 \n')
                        correct = False

                elif attribute == 'hessian_update':
                    if not str.lower(value) in ['exact','bofill']:
                        print (f'Wrong hessian_update (={value}) is given! Check the option file !!!')
                        print ('Only bofill and exact are possible options !!! Default method is bofill\n')
                        correct = False
                    else:
                        args.hessian_update = value

                elif attribute == 'reoptimize':
                    try:
                        value = int(value)
                        if not value in [0,1]:
                            print (f'Wrong reoptimize (={value}) is given! Check the option file !!!')
                            print ('Only 0 or 1 are possible. If the value is zero, given geometry is directly undergone MCD, otherwise, the molecule is reoptimized !!! Default value is 1\n')
                            correct = False
                        else:
                            args.reoptimize = value
                    except:
                        print (f'Wrong reoptimize (={value}) is given! Check the option file !!!')
                        print ('Only 0 or 1 are possible. If the value is zero, given geometry is directly undergone MCD, otherwise, the molecule is reoptimized !!! Default value is 1\n')
                        correct = False
                else:
                    wrong_attributes.append(attribute)
                    correct = False

            if len(wrong_attributes) > 0:
                content = ','.join(wrong_attributes)
                print (f'Wrong attribute(s) (={content}) is given! Check the option file !!!')
                print ('Possible attributes are \'working_directory\', \'num_relaxation\',\'step_size\',\'unit\',\'calculator\',\'command\',\'use_hessian\',\'hessian_update\',\'reoptimize\'')


    else:
        print ('option directory is not found! Default parameters are used!!!')

    return correct

def get_calculator(args):
    calculator_name = args.calculator.lower()
    if calculator_name == 'gaussian':
        from pyMCD.Calculator import gaussian
        calculator = gaussian.Gaussian(args.command)
    elif calculator_name == 'orca':
        from pyMCD.Calculator import orca
        calculator = orca.Orca()
    else:
        print (f'Wrong calculator (={calculator_name}) is given! Check the option file !!!')
        calculator = None
        return calculator
    calculator.load_content(os.path.join(args.input_directory,'qc_input'))
    basis_file = os.path.join(args.input_directory,'basis') # For Effective Core Potential
    if os.path.exists(basis_file):
        calculator.load_basis(basis_file)
    return calculator


def generate_path():
    import datetime
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_directory','-id',type=str,help='directory of inputs')
    parser.add_argument('--output_directory','-od',type=str,help='directory for saving outputs',default=None)
    parser.add_argument('--working_directory','-wd',type=str,help='scratch directory of QC programs',default='')
    parser.add_argument('--num_relaxation',type=int,help='Num relaxation for input',default=5)
    parser.add_argument('--step_size',type=float,help='Maxmial displacement',default=0.0)
    parser.add_argument('--calculator','-c',type=str,help='Name of Quantum Calculation software',default='gaussian')
    parser.add_argument('--unit','-u',type=str,help='unit',default='Hartree')
    parser.add_argument('--command',type=str,help='command for running qc package',default='g09')
    parser.add_argument('--use_hessian',type=int,default=1)
    parser.add_argument('--restart',type=int,default=0)
    parser.add_argument('--hessian_update',type=str,default='bofill')
    parser.add_argument('--reoptimize',type=int,default=1)

    args = parser.parse_args()

    use_hessian = args.use_hessian
    reoptimize = args.reoptimize

    restart = args.restart

    if use_hessian == 0:
        use_hessian = False
    else:
        use_hessian = True

    if reoptimize == 0:
        reoptimize = False
    else:
        reoptimize = True

    if restart == 0:
        restart = False
    else:
        restart = True

    input_directory = args.input_directory
    output_directory = args.output_directory
    if output_directory is None:
        output_directory = input_directory
    # If problem with input, output directory, automatically exit
    if not os.path.exists(input_directory):
        print ('Cannot find the input directory !!!')
        exit()
    elif not os.path.exists(output_directory):
        print ('Given output directory is not found !!!')
        exit()

    print (f'\ninput directory: {input_directory}')
    print (f'output directory: {output_directory}\n')

    reactant = read_reactant(input_directory) # Read geometry of reactant
    constraints, num_steps = read_bond_info(input_directory) # bond info
    correct = change_option(args) # Read option file and change values in args
    if not correct:
        exit()
    calculator = get_calculator(args) # Make calculator, you can use your own calculator!
    if calculator is None:
        exit()
    if use_hessian:
        scanner = mcd.MCD(num_relaxation = args.num_relaxation,calculator=calculator)
        scanner.hessian_update = args.hessian_update
    else:
        scanner = mcd.MCDModified(num_relaxation = args.num_relaxation,calculator=calculator)
    scanner.step_size = args.step_size
    scanner.log_directory = output_directory
    working_directory = args.working_directory
    if working_directory != '':
        if os.path.exists(working_directory):
            scanner.change_working_directory(working_directory)
        else:
            print ('working directory does not exist!!! Using output directory as default ...')
            working_directory = output_directory
    else:
        #print ('adfasdfa',output_directory)
        working_directory = output_directory
    scanner.change_working_directory(working_directory)
    print (f'working directory: {working_directory}\n')

    scanner.change_energy_unit(args.unit)
    #exit()
    # Also, write constraints information (TODO)

    num_scan = sum(num_steps.values())

    print ('\n=======================================================')
    print ('================= PYMCD RUNNING !!! ===================')
    print ('=======================================================\n')

    if use_hessian:
        pathway = scanner.scan(reactant,constraints,num_scan=num_scan,chg = reactant.chg, multiplicity = reactant.multiplicity)
    else:
        pathway = scanner.scan(reactant,constraints,num_steps,chg = reactant.chg, multiplicity = reactant.multiplicity)
    return pathway

if __name__ == '__main__':
    pathway = generate_path()     
    if len(pathway) > 0:
        print ('MCD well terminated ...')
    else:
        print ('MCD did not run properly !!!')
