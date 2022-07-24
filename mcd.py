
### python provided modules ###
import time
import os
import sys
from copy import deepcopy
import subprocess
import pickle
import datetime
import argparse

### extra common libraries ###
import numpy as np

### ace-reaction libraries ###
from pyMCD import chem
from pyMCD.utils import process
from pyMCD.utils import ic

class MCD: # MultiCoordinate Driving method for finding MEP
    
    def __init__(self,num_relaxation = 1,calculator=None):
        self.num_relaxation = num_relaxation
        self.direction = None
        self.calculator = calculator
        self.log_directory = ''
        self.energy_unit = 'Hartree'
        self.num_force_calls = 0
        self.scale = 1.0
        self.consider_redundant = False

    def __str__(self):
        content = ''
        content = content + f'num relaxation: {self.num_relaxation}\n'
        #content = content + f'working_directory: {self.working_directory}\n'
        #content = content + f'qc inputs\n\n{self.content}\n'
        return content

    def change_energy_unit(self,energy_unit):
        self.energy_unit = energy_unit
        if self.calculator is not None:
            self.calculator.energy_unit = energy_unit

    def change_working_directory(self,working_directory):
        if self.calculator is not None:
            self.calculator.change_working_directory(working_directory)
        else:
            print ('Calculator does not exist!!! Define a calcualtor first!')

    def write_log(self,content,mode='a'):
        if self.log_directory == '':
            print ('No log directory!!! We write nothing!!!')
        else:
            log_directory = os.path.join(self.log_directory,'output.log')
            with open(log_directory,mode) as f:
                f.write(content)
                f.flush()

    def update_pathway(self,calculated_molecule,index,mode='a'): # Calculated molecule that will be added to the current pathway
        save_directory = self.log_directory
        if save_directory == '':
            print ('No log directory!!! We do not save the pathway!!!')
        else:
            # First, update energy log
            with open(os.path.join(save_directory,'energy.log'),mode) as f:
                if mode == 'w':
                    f.write(f'Energy ({self.energy_unit})\n')
                f.write(f'{index} {calculated_molecule.energy}\n')
                f.flush()
            # Save xyz file
            with open(os.path.join(save_directory,'pathway.xyz'),mode) as f:
                f.write(str(index)+'\n'+str(calculated_molecule.energy)+'\n')
                for atom in calculated_molecule.atom_list:
                    content = atom.get_content()
                    f.write(content)
                f.write('\n')
                f.flush()
            # Save pickle file
            with open(os.path.join(save_directory,'pathway.pkl'),mode+'b') as f:
                data = calculated_molecule.get_minimal_data()
                data['energy'] = calculated_molecule.energy
                pickle.dump(data,f)
                f.flush()                

    def write_profile(self,trajectory):
        # Save trajectory with log file and trajectory files (xyz, pkl)
        save_directory = self.log_directory
        n = len(trajectory)
        energy_list = []
        for i in range(n):
            molecule = trajectory[i]
            energy_list.append(molecule.energy)
        n = len(energy_list)
        energy_log = open(os.path.join(save_directory,'profile.log'),'w')
        energy_log.write(f'Original Energy ({self.energy_unit}) \t Relative Energy ({self.energy_unit})\n') 
        reference_energy = energy_list[0]
        maximum_index = 0
        maximum_energy = -100000000
        for i in range(n):
            relative_energy = energy_list[i] - reference_energy
            energy = energy_list[i]
            energy_log.write(f'{energy} \t {relative_energy}\n')
            if maximum_energy < energy_list[i]:
                maximum_index = i
                maximum_energy = energy_list[i]
        # Write maxima point
        energy_log.write(f'Maxima point: {maximum_index}')
        energy_log.close()
        # Save maximal point as ts.xyz
        trajectory[maximum_index].write_geometry(os.path.join(save_directory,'ts.xyz'))
        if maximum_index == 0 or maximum_index == n - 1:
            self.write_log('Caution: Reactant/Product are found to be the highest point!\n')
         
    def get_delta_q(self,coordinate_list,constraints,num_steps):
        delta_q = dict()
        update_q = dict()
        internal_coordinates = list(constraints.keys())
        # May consider better reaction coordinate representation later!
        for constraint in constraints:
            if num_steps[constraint] == 0: # Only consider scanning parts!
                continue            
            delta_q[constraint] = constraints[constraint] - ic.get_single_q_element(coordinate_list,constraint)
            update_q[constraint] = 0
        return delta_q,update_q

    def select_coordinate(self,coordinate_list,force,delta_q,scanning_coordinates):
        # Get inner products (coefficients)
        internal_coordinates = list(delta_q.keys())
        B = ic.get_wilsonB_matrix(coordinate_list,internal_coordinates)
        B_inverse = ic.get_general_inverse(B)
        n = len(force)
        reshaped_force = -np.reshape(force,(3*n))
        dV_dq = reshaped_force@B_inverse
        minimum = 1000000
        #print (dV_dq)
        delta_es = dict()
        for idx,internal_coordinate in enumerate(scanning_coordinates):
            sign = 1
            if delta_q[internal_coordinate] < 0:
                sign = -1
            delta_e = dV_dq[idx] * sign
            delta_es[internal_coordinate] = delta_e
            if delta_e < minimum:
                minimum = delta_e
                selected_coordinate = internal_coordinate
        print ('product: ',delta_es,selected_coordinate)
        return selected_coordinate

    # Only different part between auto_ts3.py
    def update_geometry(self,molecule,constraints,num_steps,chg=None,multiplicity=None):
        num_relaxation = self.num_relaxation # N
        coordinate_list = molecule.get_coordinate_list()
        delta_q,update_q = self.get_delta_q(coordinate_list,constraints,num_steps) # Only scanning coordinates are survived
        scanning_coordinates = list(delta_q.keys())
        # If you want to put additional redundnat coordinate, activate consider_redundant, and implement your own redundant coordinate (Default: If activated, interatomic distance between bonded atoms are added)
        if self.consider_redundant:
            # Find distance that is close to each other
            adj_matrix = process.get_adj_matrix_from_distance(molecule,1.3)
            bonds = np.stack(np.where(adj_matrix>0),axis=1)
            for bond in bonds:
                tuple_bond = tuple(bond)
                if bond[0] < bond[1]:
                    if tuple_bond not in delta_q:
                        delta_q[tuple_bond] = 0.0
                    update_q[tuple_bond] = 0.0
        force = self.calculator.get_force(molecule,chg,multiplicity)
        self.num_force_calls += 1
        atom_list = molecule.atom_list
        # Make update vector
        n = len(atom_list)
        selected_coordinate = self.select_coordinate(coordinate_list,force,delta_q,scanning_coordinates)
        # Update geometry using ic module
        update_q[selected_coordinate] = delta_q[selected_coordinate]/num_steps[selected_coordinate]
        displacement = abs(update_q[selected_coordinate])
        #molecule.print_coordinate_list()
        ic.update_geometry(coordinate_list,update_q)
        num_steps[selected_coordinate] -= 1
        process.locate_molecule(molecule,coordinate_list,False)
        self.direction = selected_coordinate
        return displacement,force
        
    def get_corrected_indices(self,constraint): # In general, we use index = 1 as starting point!
        new_constraint = []
        for idx in constraint:
            new_constraint.append(idx + 1)
        new_constraint = tuple(new_constraint)
        return new_constraint

        
    def print_status(self,molecule,constraints,unit = 'degree'):
        print ('########## Current geometry info ############')
        print ('\t Target \t Current')
        coordinate_list = molecule.get_coordinate_list()
        for constraint in constraints:
            new_constraint = self.get_corrected_indices(constraint)
            target_value = constraints[constraint]
            if len(constraint) == 2:
                #print ('New',constraint,molecule.atom_list[constraint[0]].get_coordinate(),molecule.atom_list[constraint[1]].get_coordinate())
                print (f'{new_constraint}   {target_value}  {molecule.get_distance_between_atoms(constraint[0],constraint[1])}')
            elif len(constraint) == 3:
                if unit == 'degree':
                    target_value *= 180/np.pi
                print (f'{new_constraint}   {target_value}  {molecule.get_angle_between_atoms(constraint[0],constraint[1],constraint[2],unit)}')
            elif len(constraint) == 4:
                if unit == 'degree':
                    target_value *= 180/np.pi
                print (f'{new_constraint}   {target_value}  {molecule.get_dihedral_angle_between_atoms(constraint[0],constraint[1],constraint[2],constraint[3],unit)}')


    def scan(self,molecule,constraints,num_steps,chg=None,multiplicity=None): # constraints: given as bond: target_distance
        # Save zeroth point
        if type(num_steps) is int:
            new_num_steps = dict()
            for constraint in constraints:
                new_num_steps[constraint] = num_steps
            num_steps = new_num_steps # Make cubic scan coordinates
        else:
            num_steps = num_steps.copy() # It should be dict
        constraints = constraints.copy()
        ###### Write basic informations #####
        self.write_log('##### Scanning information ######\n'+str(self.calculator),'w') # Calculator information
        self.write_log(f'Num relaxation: {self.num_relaxation}\n\n')
        coordinate_list = molecule.get_coordinate_list()
        delta_q,update_q = self.get_delta_q(coordinate_list,constraints,num_steps)
        # Write reactant info
        self.write_log(f'###### Reactant information ######\n')
        if chg is not None and multiplicity is not None:
            self.write_log(f'charge: {chg}    multiplicity: {multiplicity}\n')
        else:
            self.write_log('Default charge and multiplicity values are used !!!\n')
        self.write_log(str(len(molecule.atom_list))+'\n')
        for atom in molecule.atom_list:
            self.write_log(atom.get_content())
        self.write_log('\n')
        self.write_log(f'########## Scanning coordinates ##########\n')
        for coordinate in delta_q:
            current_value = float(format(ic.get_single_q_element(coordinate_list,coordinate),'.4f'))
            target_value = constraints[coordinate]
            new_constraint = self.get_corrected_indices(coordinate)
            if len(coordinate) > 2:
                current_value *= 180/np.pi
                target_value *= 180/np.pi
            delta = round(delta_q[coordinate]/num_steps[coordinate],4)
            self.write_log(f'{new_constraint}: {current_value} -> {target_value}, {delta} angstrom per step\n')

        hostname = subprocess.check_output(['hostname'],universal_newlines=True) 
        self.write_log(f'Computing node: {hostname}\n')
        st = datetime.datetime.now()
        self.write_log(f'Starting time: {st}\n')
        print ('scanning ...')
        list_of_constraints = list(constraints.keys())
        self.write_log('Set up done!!!\n')
        self.write_log('Initialization, Relaxing ....\n')
        calculated_data = self.calculator.relax_geometry(molecule,list_of_constraints,chg,multiplicity,'test',None,1000,True)
        try:
            self.num_force_calls += len(calculated_data.atomcoords)
        except:
            a = 1
        self.write_log('Relaxation finished! Start scanning ....\n')
        #energy = self.get_energy(molecule,chg,multiplicity)
        copied_molecule = molecule.copy()
        copied_molecule.energy = molecule.energy
        trajectory = [copied_molecule] # Trajectory is list of coordinates
        #self.print_status(molecule,constraints)
        self.direction = list(constraints.keys())[0] # Just default direction, it's changed in update_geometry
        energy_list = [molecule.energy]
        self.update_pathway(copied_molecule,0,'w')
        total_num_scans = sum(num_steps.values())
        #force_list = []
        original_num_steps = num_steps.copy()
        digit = 3
        while sum(num_steps.values()) > 0:
            starttime = datetime.datetime.now()
            #self.write_log(f'[{starttime}] Performing one step search ... ({len(trajectory)-1}/{total_num_scans}) completed ... \n')
            # Update/Relax geometry
            self.print_status(molecule,constraints)
            list_of_constraints = []
            for constraint in num_steps:
                if num_steps[constraint] > 0:
                    list_of_constraints.append(constraint)
            displacement,force = self.update_geometry(molecule,constraints,num_steps,chg,multiplicity)
            num_force_calls = 1
            selected_coordinate = self.direction
            tmp_st = str(datetime.datetime.now())[:-digit]
            content = f'[{tmp_st}] Progress: '
            for coordinate in original_num_steps:
                total_num = original_num_steps[coordinate]
                progress_num = total_num - num_steps[coordinate]
                new_constraint = self.get_corrected_indices(coordinate)
                content = content + f' {new_constraint}: {progress_num}/{total_num}'
            #content = content + f' Selected coordinate: {selected_coordinate}'
            self.write_log(content+'\n')
            #force_list.append(force)
            print ('status',displacement,num_steps)
            calculated_data = self.calculator.relax_geometry(molecule,list_of_constraints,chg,multiplicity,'test',self.num_relaxation,displacement*self.scale,True)
            scfenergies = calculated_data.scfenergies
            print ('relaxation energy: ',scfenergies[0] - molecule.energy)
            try:
                num_force_calls += len(calculated_data.atomcoords)
            except:
                a = 1
            endtime = datetime.datetime.now()
            energy = molecule.energy
            copied_molecule = molecule.copy()
            copied_molecule.energy = energy
            # Summary log for a single step
            delta_e = energy - energy_list[-1]
            delta_time = endtime - starttime
            word = 'Increased'
            exponent = int(np.log10(np.abs(delta_e)))
            x = delta_e * 10**(-exponent)
            if x < 0:
                word = 'Decreased'           
                x *= -1 
            x = format(x,'.4f')
            self.num_force_calls += num_force_calls
            endtime = str(endtime)[:-digit]
            delta_time = str(delta_time)[:-digit]
            self.write_log(f'[{endtime}] {x}E{exponent} {self.energy_unit} has {word} after {num_force_calls} force calls! {delta_time} Taken ...\n')
            energy_list.append(energy)
            # Save/Copy new molecule
            self.update_pathway(copied_molecule,len(trajectory),'a')
            trajectory.append(copied_molecule)
            delta_energy = energy_list[-1] - energy_list[-2]
            print ('Energy increase: ',delta_energy/displacement)
            print ('Energy: ',energy_list)
        self.write_log(f'[{datetime.datetime.now()}] Scan completed ...\n')
        self.write_log(f'Total {self.num_force_calls} force calls performed ...\n')
        self.direction = None
        self.write_profile(trajectory)
        et = datetime.datetime.now()
        self.write_log(f'End time: {et}\n')
        self.write_log(f'Taken time: {et-st}\n')
        self.calculator.clean_scratch()
        return trajectory

