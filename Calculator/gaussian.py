### python provided modules ###
import time
import os
import sys
from copy import deepcopy
import subprocess
import pickle
import datetime
import argparse
import numpy as np


### Module for reading gaussian files ###
import cclib

### ace-reaction libraries ###
from pyMCD.utils import process


class Gaussian:

    def __init__(self,command='g09'):
        self.working_directory = os.getcwd()
        self.command = command
        self.content='#N pm6 scf(xqc) '
        self.energy_unit = 'Hartree'
        self.basis = ''

    def __str__(self):
        content = ''
        content = content + f'working_directory: {self.working_directory}\n'
        content = content + f'command: {self.command}\n'
        content = content + f'Energy: {self.energy_unit}\n\n'
        content = content + f'###### qc_input ######\n{self.content}\n'
        return content

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

    def change_working_directory(self,working_directory):
        # Get current reaction coordinate
        if not os.path.exists(working_directory):
            print ('Working directory does not exist! Creating new directory ...')
            try:
                os.system(f'mkdir {self.working_directory}')
            except:
                print ('Cannot create working directory!!!\n Recheck your working directory ...')
                working_directory = self.working_directory
        self.working_directory = working_directory
        os.environ['GAUSS_SCRDIR'] = working_directory

    def get_content(self):
        return self.content 

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

    def write_molecule_info(self,f,molecule,chg,multiplicity):
        f.write(f'{chg} {multiplicity}\n')
        for atom in molecule.atom_list:
            f.write(atom.get_content())
        f.write('\n')
    
    def write_basis(self,f):
        if self.basis != '':
            f.write(f'{self.basis}\n')

    def make_input(self,molecule,chg=None,multiplicity=None,option='sp',file_name='test',constraints={},extra=''):
        f = open(f'{file_name}.com','w')
        if chg is None or multiplicity is None:
            chg, multiplicity = self.get_default_mol_params(molecule)
        content = self.get_content()
        content = content + option + f'{extra}\n\ntest\n\n'
        f.write(content)
        self.write_molecule_info(f,molecule,chg,multiplicity)
        ### Write constraints
        for constraint in constraints:
            constraint_info = [] 
            if len(constraint) == 2:
                constraint_info.append('B')
            elif len(constraint) == 3:
                constraint_info.append('A')
            else:
                constraint_info.append('D') 
            for index in constraint:
                constraint_info.append(str(index+1))
            constraint_info.append('F')
            f.write(' '.join(constraint_info)+'\n')
        f.write('\n') # Additional empty line required
        self.write_basis(f)
        f.write('\n')
        f.close()


    def get_energy(self,molecule,chg=None,multiplicity=None,file_name='test',extra=''):
        current_directory = os.getcwd()
        os.chdir(self.working_directory)
        self.make_input(molecule,chg,multiplicity,option='sp',file_name=file_name,extra=extra)
        os.system(f'{self.command} {file_name}.com')
        # Read output
        p = cclib.parser.Gaussian('{file_name}.log')
        data = p.parse()
        converter = 1
        os.chdir(current_directory)
        if self.energy_unit == 'kcal':
            converter = 23.06
        if self.energy_unit == 'Hartree':
            converter = 0.036749326681
        #os.system('mv new.chk old.chk')
        os.chdir(current_directory)
        return converter*data.scfenergies[-1]

    def get_force(self,molecule,chg=None,multiplicity=None,file_name='test',extra=' Symmetry = None'):
        current_directory = os.getcwd()
        os.chdir(self.working_directory)
        self.make_input(molecule,chg,multiplicity,option='force',file_name=file_name,extra=extra)
        os.system(f'{self.command} {file_name}.com')
        # Read output
        p = cclib.parser.Gaussian(f'{file_name}.log')
        data = p.parse()
        #os.system('mv new.chk old.chk')
        os.chdir(current_directory)
        return data.grads[-1]

    def optimize_geometry(self,molecule,constraints={},chg=None,multiplicity=None,file_name='test',extra='',return_data = False):
        current_directory = os.getcwd()
        os.chdir(self.working_directory)
        self.make_input(molecule,chg,multiplicity,'opt',file_name,constraints,extra)
        os.system(f'{self.command} {file_name}.com')
        # Read output
        p = cclib.parser.Gaussian(f'{file_name}.log')
        data = p.parse()
        # Get minimal energy geometry
        index = np.argmin(data.scfenergies)
        if index > len(data.atomcoords) - 1:
            index = -1
        coordinate_list = data.atomcoords[index] # Sometimes length is not matched, it performs extra scfenergy calculation
        converter = 1
        if self.energy_unit == 'kcal':
            converter = 23.06
        if self.energy_unit == 'Hartree':
            converter = 0.036749326681
            #converter = 1/27.2114
        energy = data.scfenergies[index] * converter
        process.locate_molecule(molecule,coordinate_list,False)
        molecule.energy = energy
        #os.system('mv new.chk old.chk')
        os.chdir(current_directory)
        if return_data:
            return data
                

class GaussianMCD(Gaussian):
    
    def __init__(self,command='g09'):
        super().__init__(command)
        self.old_chk_name = 'old'
        self.new_chk_name = 'new'

    
    def relax_geometry(self,molecule,constraints,chg=None,multiplicity=None,file_name='test',num_relaxation=5,maximal_displacement=1000,return_data=False):
        if maximal_displacement < 100:
            max_step = int(maximal_displacement*100) + 1
            extra = f'(modredundant,loose,maxcycles={num_relaxation},maxstep={max_step},notrust) Symmetry=None'
        else:
            extra = f'(modredundant,loose,maxcycles=15) Symmetry = None'
        data = self.optimize_geometry(molecule,constraints,chg,multiplicity,file_name,extra,return_data)
        if return_data:
            return data

    def get_content(self):
        # Check whether old one exists
        content = self.content 
        if os.path.exists(f'{self.old_chk_name}.chk'):
            print ('check point exists!!!')
            content = f'%oldchk={self.old_chk_name}.chk\n' + content
            content = content + ' guess=read '
        content = f'%chk={self.new_chk_name}.chk\n' + content        
        return content      

    def change_chk_file(self):
        old_chk_directory = os.path.join(self.working_directory,self.old_chk_name+'.chk')
        new_chk_directory = os.path.join(self.working_directory,self.new_chk_name+'.chk')
        os.system(f'mv {new_chk_directory} {old_chk_directory}')

    def get_energy(self,molecule,chg=None,multiplicity=None,file_name='test',extra=''):
        energy = super().get_energy(molecule,chg,multiplicity,file_name,extra)
        self.change_chk_file()
        return energy

    def get_force(self,molecule,chg=None,multiplicity=None,file_name='test',extra=' Symmetry = None'):
        force = super().get_force(molecule,chg,multiplicity,file_name,extra)
        self.change_chk_file()
        return force

    def optimize_geometry(self,molecule,constraints = {},chg=None,multiplicity=None,file_name='test',extra='',return_data = False):
        data = super().optimize_geometry(molecule,constraints,chg,multiplicity,file_name,extra,return_data)
        self.change_chk_file()
        if return_data:
            return data

    def clean_scratch(self,file_name='test.com'):
        working_directory = self.working_directory
        chk_directory = os.path.join(working_directory,'old.chk')
        os.system(f'rm {chk_directory}')
        chk_directory = os.path.join(working_directory,'new.chk')
        os.system(f'rm {chk_directory}')
        file_directory = os.path.join(working_directory,'test.com')
        os.system(f'rm {file_directory}')
        file_directory = os.path.join(working_directory,'test.log')
        os.system(f'rm {file_directory}')
