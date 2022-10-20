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

### ace-reaction libraries ###
from pyMCD.utils import process

'''
You can define your own calculator, depending on your using software!!!
Here, you must note that below functions must be defined in order to run calculations with MCD: __init__, get_energy, get_force, relax_geometry
See orca.py and gaussian.py as examples. How they
You can also try making calculations using ase. In our case, we did not use ASE modules due to some optimizing issues ...
'''
class Calculator:

    def __init__(self):
        self.working_directory = os.getcwd()
        self.energy_unit = 'Hartree'

    def change_working_directory(self,working_directory):
        # Get current reaction coordinate
        if not os.path.exists(working_directory):
            print ('Working directory does not exist! Creating new directory ...')
            os.system(f'mkdir {self.working_directory}')
            self.working_directory = working_directory

    def load_content(self,template_directory):
        '''
        Load route section

        For examples, self.contents will be "#N b3lyp/6-31g(d)" for gaussian calculation
        '''
        content = ''
        with open(template_directory) as f:
            for line in f:
                content = content + line
        self.content = content

    def load_basis(self,basis_directory):
        """
        Load additional basis information(Effective core ~~)
        """
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


    def get_energy(self,molecule,chg=None,multiplicity = None):
        '''
        Must return energy with desired unit defined in the Calculator
        '''
        pass


    def get_force(self,molecule,chg=None,multiplicity=None):
        '''
        Must return force with desired unit defined in the Calculator
        '''
        pass


    def relax_geometry(self,molecule,constraints,chg=None,multiplicity=None,num_relaxation=5,maximal_displacement=1000,return_data=False):
        '''
        Optimize the geometry under constraints for given relaxation steps.
        '''
        pass
