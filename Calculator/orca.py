### python provided modules ###
import time
import os
import sys
from copy import deepcopy
import subprocess
import pickle
import datetime
import argparse
import distutils.spawn

import numpy as np

### ace-reaction libraries ###
from pyMCD.utils import process

'''
You can define your own calculator, depending on your using software!!!
Here, you must note that below functions must be defined in order to run calculations with MCD: __init__, get_energy, get_force, relax_geometry
See orca.py and gaussian.py as examples. How they
You can also try making calculations using ase. In our case, we did not use ASE modules due to some optimizing issues ...
'''
class Calculator:

    def __init__(self,command='orca'):
        self.working_directory = os.getcwd()
        check = distutils.spawn.find_executable(command)
        if check is None:
            print ('orca not found!')
            exit()
        self.command = command
        self.energy_unit = 'Hartree'
        self.nproc=1
        self.content = '! b3lyp 6-31g(d)'

    def change_working_directory(self,working_directory):
        # Get current reaction coordinate
        if not os.path.exists(working_directory):
            print ('Working directory does not exist! Creating new directory ...')
            os.system(f'mkdir {self.working_directory}')
            self.working_directory = working_directory

    def load_content(self,template_directory):
        content = ''
        with open(template_directory) as f:
            for line in f:
                content = content + line
        self.content = content

    def load_basis(self,basis_directory):
        basis = ''
        with open(basis_directory) as f:
            for line in f:
                basis = basis + line
        self.basis = basis

    def get_default_mol_params(self,molecule):
        try:
            chg = molecule.get_chg()
        except:
            chg = 0
        try:
            e_list = molecule.get_num_of_lone_pair_list()
            num_of_unpaired_e = len(np.where((2*e_list) % 2 == 1)[0])
            multiplicity = num_of_unpaired_e + 1
        except:
            z_sum = np.sum(molecule.get_z_list())
            multiplicity = (z_sum - chg) % 2 + 1
        return chg,multiplicity

    def write_input_file(self,molecule,chg, multiplicity, constraints={},option ='sp',extra=''):

        if chg is None or multiplicity is None:
            chg, multiplicity = self.get_default_mol_params(molecule)

        inpstring = self.content.strip()
        if option == 'force':
            inpstring += ' EnGrad\n\n'  # SOSCF SlowConv \n\n'
        elif option == 'sp':
            inpstring += '\n\n'
        elif option == 'opt':
            inpstring += ' OPT\n\n'  # SOSCF SlowConv \n\n'
        '''
        inpstring += '%scf\nMaxIter 300\nconvergence strong\n sthresh 1e-7\n'
        inpstring += 'thresh 1e-11\n tcut 1e-13 \n directresetfreq 1 \n SOSCFStart 0.00033\nend\n'
        inpstring += '%scf\nMaxIter 300\nend\n'
        inpstring += '\n%maxcore 1000\n\n'
        inpstring += '%pal\nnproc {}\nend\n\n'.format(self.nproc)
        '''
        inpstring += '%geom\nConstraints\n'
        for constraint in constraints:
            constraint_info = []
            if len(constraint) == 2:
                constraint_info.append('{B')
            elif len(constraint) == 3:
                constraint_info.append('{A')
            else:
                constraint_info.append('{D')
            for index in constraint:
                constraint_info.append(str(index-1))
            constraint_info.append('C}')
            inpstring +=  (' '.join(constraint_info)+'\n')
        inpstring += 'end\n'
        inpstring += extra
        inpstring += 'end\n'

        inpstring += '\n*xyz {} {}\n'.format(chg, multiplicity)
        for atom in molecule.atom_list:
            inpstring += atom.get_content()

        inpstring += '*'
        tempfilename = 'tempORCAinp_{}'.format(multiplicity)
        tempfile = open(tempfilename, 'w')
        tempfile.write(inpstring)
        tempfile.close()
        return tempfilename

    def get_energy(self,molecule,chg=None,multiplicity = None):
        '''
        Must return energy with desired unit defined in the Calculator
        '''
        # Write input file
        tempfilename = self.write_input_file(molecule, chg, multiplicity,option='sp')
        runscr = self.working_directory
        os.system('mkdir -p {}'.format(runscr))
        os.system('mv {} {}/'.format(tempfilename, runscr))
        cmd = 'cd {}; {} {} > {}/{}.log'.format(runscr, self.command, tempfilename, runscr, tempfilename)
        os.system(cmd)

        energy =self.parse_energy(multiplicity, runscr, tempfilename)

        return energy


    def get_force(self,molecule,chg=None,multiplicity=None):
        '''
        Must return force with desired unit defined in the Calculator
        '''
        # Write input file
        tempfilename = self.write_input_file(molecule, chg, multiplicity,option='force')
        runscr = self.working_directory
        os.system('mkdir -p {}'.format(runscr))
        os.system('mv {} {}/'.format(tempfilename, runscr))
        cmd = 'cd {}; {} {} > {}/{}.log'.format(runscr, self.command, tempfilename, runscr, tempfilename)
        os.system(cmd)

        force =self.parse_force(multiplicity, runscr, tempfilename)

        return -force

    def parse_energy(self, multiplicity,runscr,tempfilename):
        engradpath = runscr+'/{}.energy'.format(tempfilename)
        with open(engradpath) as engradfile:
            engradlines = engradfile.readlines()

        temp = 100000
        for i, lines in enumerate(engradlines):
            if '$energy\n' in lines:
                temp = i
            if i > temp+1:
                energy = float(lines.split()[-1])
                break
        return energy

    def parse_force(self, multiplicity,runscr,tempfilename):
        engradpath = runscr+'/{}.engrad'.format(tempfilename)
        with open(engradpath) as engradfile:
            engradlines = engradfile.readlines()

        temp = 100000
        tmp = []
        tmp2 = []
        for i, lines in enumerate(engradlines):
            if '# The current gradient in Eh/bohr\n' in lines:
                temp = i
            if i > temp+1:
                if "#" in lines:
                    break
                tmp2.append(float(lines.split()[0]))
            if len(tmp2) == 3:
                tmp.append(tmp2)
                tmp2 = []
        return np.asarray(tmp)

    def parse_opt(self, multiplicity,runscr,tempfilename):
        engradpath = runscr+'/{}.xyz'.format(tempfilename)
        with open(engradpath) as engradfile:
            engradlines = engradfile.readlines()
        xyz = []
        for i, lines in enumerate(engradlines[2:]):
            a,x,y,z = lines.strip().split()
            xyz.append([x,y,z])
        xyz = np.array(xyz,dtype='float32')

        engradpath = runscr+'/{}.energy'.format(tempfilename)
        with open(engradpath) as engradfile:
            engradlines = engradfile.readlines()

        temp = 100000
        for i, lines in enumerate(engradlines):
            if '$energy\n' in lines:
                temp = i
            if i > temp:
                print(lines)
                energy = float(lines.split()[-1])
                break

        return xyz, energy


    def optimize_geometry(self,molecule,constraints={},chg=None,multiplicity=None,extra='',return_data = False):
        current_directory = os.getcwd()
        # Write input file
        tempfilename = self.write_input_file(molecule, chg, multiplicity,constraints,option='opt',extra=extra)
        runscr = self.working_directory
        os.system('mkdir -p {}'.format(runscr))
        os.system('mv {} {}/'.format(tempfilename, runscr))
        cmd = 'cd {}; {} {} > {}/{}.log'.format(runscr, self.command, tempfilename, runscr, tempfilename)
        os.system(cmd)
        xyz, energy = self.parse_opt(multiplicity,runscr,tempfilename)
        process.locate_molecule(molecule,xyz,False)
        molecule.energy = energy
        os.chdir(current_directory)
        if return_data:
            return xyz

    def relax_geometry(self,molecule,constraints,chg=None,multiplicity=None,num_relaxation=5,maximal_displacement=1000,return_data=False):
        '''
        '''
        maximal_displacement=-0.3
        if num_relaxation is None:
            extra = f'Trust {maximal_displacement}\n'
        else:
            extra = f'Trust {maximal_displacement}\nMaxIter {num_relaxation}\n'
        data = self.optimize_geometry(molecule,constraints,chg,multiplicity,extra,return_data)
        return data
