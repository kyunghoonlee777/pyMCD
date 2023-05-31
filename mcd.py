
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
import scipy
from scipy.optimize import minimize
from scipy.optimize import Bounds
from scipy.signal import savgol_filter 


### ace-reaction libraries ###
from pyMCD import chem
from pyMCD.utils import process
from pyMCD.utils import ic

class MCD: # MultiCoordinate Driving method for finding MEP
    
    def __init__(self,num_relaxation = 3,calculator=None):
        self.num_relaxation = num_relaxation
        self.calculator = calculator
        self.log_directory = ''
        self.energy_unit = 'Hartree'
        self.num_force_calls = 0
        self.num_hessian_calls = 0
        self.step_size = 0.0
        self.use_hessian = False
        self.hessian_update = 'Bofill'
        self.num_filter_iteration = 10
        self.previous_gradient = None
        self.previous_hessian = None
        self.previous_geometry = None


    def __str__(self):
        content = ''
        content = content + f'num relaxation: {self.num_relaxation}\n'
        content = content + f'step size: {self.step_size}\n'
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

    def update_pathway(self,calculated_molecule,mode='a'): # Calculated molecule that will be added to the current pathway
        save_directory = self.log_directory
        if save_directory == '':
            print ('No log directory!!! We do not save the pathway!!!')
        else:
            # Save xyz file
            with open(os.path.join(save_directory,'pathway.xyz'),mode) as f:
                f.write(str(len(calculated_molecule.atom_list))+'\n'+str(calculated_molecule.energy)+'\n')
                for atom in calculated_molecule.atom_list:
                    content = atom.get_content()
                    f.write(content)
                f.write('\n')
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
        
        # Always make it into kcal/mol
        converter = 1
        if self.energy_unit == 'Hartree':
            converter = 627.5094740631 
        elif self.energy_unit == 'eV':
            converter = 23.0605506577
        energy_log.write(f'Original Energy ({self.energy_unit}/mol) \t Relative Energy (kcal/mol)\n') 
        reference_energy = energy_list[0]
        maximum_index = 0
        maximum_energy = -100000000
        relative_energy_list = []
        for i in range(n):
            relative_energy = energy_list[i] - reference_energy
            relative_energy_list.append(relative_energy*converter)
            energy = energy_list[i]
            energy_log.write(f'{energy} \t {relative_energy*converter}\n')
            if energy > maximum_energy:
                maximum_energy = energy
                maximum_index = i
        # Save figure
        self.save_figure(relative_energy_list,'kcal')
        # Write maxima point ...
        trajectory[maximum_index].write_geometry(os.path.join(save_directory,f'ts.xyz'))
        
        # Search valley indices with using Golay Filter
        for i in range(self.num_filter_iteration):
            energy_list = savgol_filter(energy_list,9,4) 

        left_min = None
        max_index = None
        ts_indices = []
        right_min = None
        max_energy_index = None
        max_energy = -1000000

        # Assume that curve is now well smoothed ...
        for i in range(n-1):
            #print (left_min,max_index)
            if left_min is None:
                if energy_list[i] < energy_list[i+1]: # Increasing
                    left_min = i
            else:
                if max_index is None:
                    if energy_list[i] > energy_list[i+1]: # Decreasing
                        max_index = i
                else:
                    if energy_list[i] < energy_list[i+1]: # If reincreasing ...
                        right_min = i
                        ts_indices.append(max_index)
                        left_min = i
                        max_index = None # Search for the next local maximum
                    elif i == n-2:
                        right_min = i
                        ts_indices.append(max_index)

        # Write maxima points
        if len(ts_indices) > 0:
            content = ', '.join([str(ts_index) for ts_index in ts_indices])
            energy_log.write(f'\nLocal maxima around {content}th point(s) are found !!!\n\n')
        else:
            energy_log.write(f'No local maxima points are found !!! May be barrierless reaction ???\n\n')
        energy_log.close()

        if len(ts_indices) == 0:
            self.write_log('Caution: Barrierless reaction is found !!!\n')

    def save_figure(self,data,unit='kcal'):
        try:
            import matplotlib.pyplot as plt
        except:
            return None
        label_fontsize = 12
        tick_fontsize = 10
        ax = plt.subplot(111)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        x = list(range(len(data)))
        plt.plot(x,np.array(data),'s-',markersize=4,color='black')
        plt.xlabel('Step',fontsize=label_fontsize)
        plt.ylabel(f'Energy ({unit}/mol)',fontsize=label_fontsize)
        plt.xticks(fontsize=tick_fontsize)
        plt.yticks(fontsize=tick_fontsize)
        plt.savefig(os.path.join(self.log_directory,'profile.png'))

         
    def get_delta_q(self,coordinate_list,constraints):
        delta_q = dict()
        update_q = dict()
        internal_coordinates = list(constraints.keys())
        # May consider better reaction coordinate representation later!
        for constraint in constraints:
            delta_q[constraint] = constraints[constraint] - ic.get_single_q_element(coordinate_list,constraint)
            update_q[constraint] = 0
        return delta_q,update_q


    def find_solution(self,grad,hessian,update_size,delta_q,scanning_coordinates):
        # grad: n x 1, hessian: n x n
        # Make constraints
        x0 = [] # Initial guess
        lower_bounds = []
        upper_bounds = []
        jac = []
        b = np.squeeze(grad) # M x 1 -> M
        H = hessian
        n = len(scanning_coordinates)
        constant = 0
        # Ao -> Bohr
        for constraint in scanning_coordinates:
            delta = delta_q[constraint]
            if delta > 0:
                # Upper bound
                upper_bound = min(update_size,delta)
                lower_bounds.append(0)
                upper_bounds.append(upper_bound)
                value = np.random.rand(1)[0]
                constant += value
                x0.append(value)
                jac.append(1.0)
            else:
                # Lower bound
                lower_bound = max(-update_size,delta)
                lower_bounds.append(lower_bound)
                upper_bounds.append(0)
                value = np.random.rand(1)[0]
                constant += value
                x0.append(-value)
                jac.append(-1.0)
        jac = np.array(jac)
        x0 = np.array(x0)
        #print ('cons',constant)
        x0 /= constant 
        x0 *= update_size
        #print ('x0',x0)
        constraints = [{'type':'eq','fun':lambda x: np.sum(x*jac) - update_size,'jac':lambda x:jac}]
        bounds = Bounds(lower_bounds,upper_bounds)
        energy_function = lambda x: 1/2*x.T@H@x + np.sum(x*b)
        jac_energy_function = lambda x: x@H + b
        #hessian_energy_function = lambda x: H
        res = minimize(energy_function,x0,method='SLSQP',jac=jac_energy_function,constraints=constraints,bounds=bounds)
        solution = dict()
        #print ('energy increase: ',energy_function(res.x))
        print (H, b, res.x)
        print ('Expected energy change:',energy_function(res.x))

        for i,constraint in enumerate(scanning_coordinates):
            solution[constraint] = res.x[i]
        return solution

    def select_coordinate(self,coordinate_list,dV_dx,delta_q,scanning_coordinates):
        # Get inner products (coefficients)
        constraints = list(delta_q.keys())
        dV_dq = ic.get_gradient_q(dV_dx,coordinate_list,constraints)        
        minimum = 1000000
        #print (dV_dq)
        delta_es = dict()
        for idx,internal_coordinate in enumerate(constraints):
            if internal_coordinate not in scanning_coordinates:
                continue
            sign = 1
            if delta_q[internal_coordinate] < 0:
                sign = -1
            delta_e = dV_dq[idx][0] * sign
            delta_es[internal_coordinate] = delta_e
            if delta_e < minimum:
                minimum = delta_e
                selected_coordinate = internal_coordinate
        #print ('product: ',delta_es,selected_coordinate)
        print ('Expected energy change:',minimum)
        return selected_coordinate


    def get_gradient(self,molecule,chg=None,multiplicity=None):
        n = len(molecule.atom_list)
        force = self.calculator.get_force(molecule,chg,multiplicity)
        dV_dx = -np.reshape(force,(3*n,1))
        self.num_force_calls += 1
        return dV_dx


    def get_hessian(self,molecule,chg=None,multiplicity=None):
        n = len(molecule.atom_list)
        if self.previous_hessian is None:
            force,hessian = self.calculator.get_hessian(molecule,chg,multiplicity)
            dV_dx = -np.reshape(force,(3*n,1))
            self.num_hessian_calls += 1
            return dV_dx, hessian
        else:
            hessian_update = str.lower(self.hessian_update)
            if hessian_update == 'bofill':
                dV_dx = self.get_gradient(molecule,chg,multiplicity)
                G_i = self.previous_hessian
                dg_i = dV_dx - self.previous_gradient
                x_i = np.reshape(molecule.get_coordinate_list(),(3*n,1))
                x_i_1 = np.reshape(self.previous_geometry,(3*n,1))
                dx_i = x_i - x_i_1 
                #print (np.linalg.norm(dx_i))
                E_i = dg_i - G_i@dx_i
                hessian_MS = G_i + (E_i@E_i.T)/(E_i.T@dx_i)
                hessian_PSB = G_i + (E_i@dx_i.T+dx_i@E_i.T)/(dx_i.T@dx_i) - (dx_i.T@E_i)*(dx_i@dx_i.T)/(dx_i.T@dx_i)**2
                gamma = 1 - (dx_i.T@E_i)**2/((dx_i.T@dx_i)*(E_i.T@E_i))
                hessian = (1-gamma) * hessian_MS + gamma * hessian_PSB
                return dV_dx,hessian
            elif hessian_update == 'exact':
                force,hessian = self.calculator.get_hessian(molecule,chg,multiplicity)
                dV_dx = -np.reshape(force,(3*n,1))
                self.num_hessian_calls += 1
                return dV_dx, hessian


    # Only different part between auto_ts3.py
    def update_geometry(self,molecule,constraints,updates,chg=None,multiplicity=None):
        num_relaxation = self.num_relaxation # N
        coordinate_list = molecule.get_coordinate_list()
        delta_q,update_q = self.get_delta_q(coordinate_list,constraints) # Only scanning coordinates are survived
        # If you want to put additional redundnat coordinate, activate consider_redundant, and implement your own redundant coordinate (Default: If activated, interatomic distance between bonded atoms are added)
        scanning_coordinates = []
        atom_list = molecule.atom_list
        n = len(atom_list)
        for constraint in delta_q:
            if abs(delta_q[constraint]) > 0.001:
                scanning_coordinates.append(constraint)
        if len(scanning_coordinates) > 1:
            update_size = updates[scanning_coordinates[0]]
            # Make update vector
            if self.use_hessian:
                dV_dx,hessian = self.get_hessian(molecule,chg,multiplicity)
                dV_dq = ic.get_gradient_q(dV_dx,coordinate_list,constraints)
                H_q = ic.get_hessian_q(dV_dx,hessian,coordinate_list,constraints)
                
                # Make reduce function
                reduced_indices = []
                for idx, internal_coordinate in enumerate(constraints):
                    if internal_coordinate in scanning_coordinates:
                        reduced_indices.append(idx)
                reduce_function = np.ix_(reduced_indices,reduced_indices)
                dV_dq = dV_dq[reduced_indices]
                H_q = H_q[reduce_function]
                #print (dV_dq)
                #print (hessian,H_q)
                # Make update vector
                solution = self.find_solution(dV_dq,H_q,update_size,delta_q,scanning_coordinates)
                for constraint in solution:
                    update_q[constraint] = solution[constraint]

                # After using them, update previous ...
                self.previous_gradient = dV_dx
                self.previous_hessian = hessian
                self.previous_geometry = molecule.get_coordinate_list()
            else:
                dV_dx = self.get_gradient(molecule,chg,multiplicity)
                selected_coordinate = self.select_coordinate(coordinate_list,dV_dx,delta_q,scanning_coordinates)
                # Update geometry using ic module
                update_q[selected_coordinate] = updates[selected_coordinate]
        else:
            constraint = scanning_coordinates[0]
            delta = delta_q[constraint]
            if self.use_hessian:
                update_size = updates[constraint]
                if delta > 0:
                    update_q[constraint] = min(delta,update_size)
                else:
                    update_q[constraint] = max(delta,-update_size)
            else:
                update_q[constraint] = updates[constraint]
            force = None
        #molecule.print_coordinate_list()
        ic.update_geometry(coordinate_list,update_q) 
        process.locate_molecule(molecule,coordinate_list)
        
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
            current_value = molecule.get_internal_coordinate(constraint,unit)
            current_value = float(format(current_value,'.4f'))
            if len(constraint) > 2:
                if unit == 'degree':
                    target_value *= 180/np.pi
            print (f'{new_constraint}   {target_value}  {current_value}')

    def make_restart_file(self,molecule,constraints,num_scan):
        save_directory = self.log_directory
        if save_directory == '':
            print ('No log directory!!! We do not save the pathway!!!')
        else:
            # Save xyz file
            with open(os.path.join(save_directory,'new_R.com'),'w') as f:
                f.write('{molecule.get_chg()} {molecule.get_multiplicity()}\n')
                for atom in calculated_molecule.atom_list:
                    f.write(atom.get_content())
                f.write('\n')
                f.flush()

            with open(os.path.join(save_directory,'new_coordinates'),'w') as f:
                for i,constraint in enumerate(constraints):
                    value = constraints[constraint]
                    new_constraint = [str(index+1) for index in constraint]
                    content = ' '.join(new_constraint)
                    if i == len(constraints) - 1:
                        num_step = num_scan
                    else:
                        num_step = 0
                    f.write(f'{content} {value} {num_step}\n')
      


    def scan(self,molecule,constraints,num_steps,chg=None,multiplicity=None,restart = False,reoptimize = True): # constraints: given as bond: target_distance

        if type(num_steps) is int:
            new_num_steps = dict()
            for constraint in constraints:
                new_num_steps[constraint] = num_steps
            num_steps = new_num_steps # Make cubic scan coordinates
        else:
            num_steps = num_steps.copy() # It should be dict
        num_scan = sum(num_steps.values())

        self.calculator.set_error_directory(self.log_directory)
        constraints = constraints.copy()
        coordinate_list = molecule.get_coordinate_list()
        delta_q,update_q = self.get_delta_q(coordinate_list,constraints)
        # Write reactant info
        self.write_log(f'###### Reactant information ######\n','w')
        if chg is not None and multiplicity is not None:
            self.write_log(f'charge: {chg}    multiplicity: {multiplicity}\n')
        else:
            self.write_log('Default charge and multiplicity values are used !!!\n')
        self.write_log(str(len(molecule.atom_list))+'\n')
        for atom in molecule.atom_list:
            self.write_log(atom.get_content())
        self.write_log('\n')
        self.write_log(f'######### Scanning coordinates #########\n')

        # Make updates
        updates = dict()
        maximum_delta = -100
        total_delta = 0
        digit = 3
        for coordinate in delta_q:
            current_value = ic.get_single_q_element(coordinate_list,coordinate)
            target_value = constraints[coordinate]
            new_constraint = self.get_corrected_indices(coordinate)
            if len(coordinate) > 2:
                current_value *= 180/np.pi
                target_value *= 180/np.pi
            current_value = float(format(current_value,'.4f'))
            target_value = float(format(target_value,'.4f'))
            delta = delta_q[coordinate]
            total_delta += abs(delta)
            delta /= num_steps[coordinate]
            updates[coordinate] = delta
            if abs(delta) > maximum_delta:
                maximum_delta = abs(delta)
            self.write_log(f'{coordinate}: \t {current_value} -> {target_value}\n')
        
        # Reset updates when hessian is used ...
        if self.use_hessian:
            for constraint in constraints:
                updates[constraint] = total_delta/num_scan

        # Change step size

        self.write_log('\n ################## Update Information ################### \n')
        self.write_log(f'Use hessian: {self.use_hessian}\n')
        self.write_log(f'Hessian update (if used): {self.hessian_update}\n')

        original_step_size = self.step_size
        # Set default step size: (Maximum delta x/num_relaxation * 3)
        if self.step_size == 0:
            if self.use_hessian:
                self.step_size = total_delta/num_scan
                self.step_size *= 3/self.num_relaxation
            else:
                self.step_size = maximum_delta/self.num_relaxation * 3
            self.write_log(f'Step size is set to {round(self.step_size,4)} as the default value!!!\n\n')

        # Set default step size: (Maximum delta x/num_relaxation * 3)
        self.write_log(f'########## Scanning information ##########\n'+str(self)+'\n')
        ###### Write basic informations #####
        self.write_log('##### Calculator information ######\n'+str(self.calculator)) # Calculator information
        self.write_log(f'\nAll electronic energies are written in {self.energy_unit}! \n')
        self.write_log(f'All distances are written in angstrom! \n')
        self.write_log(f'All angles are written in degree(o)! \n\n')

        st = datetime.datetime.now()
        self.write_log(f'Starting time: {st}\n')
        print ('scanning ...')
        list_of_constraints = list(constraints.keys())
        self.write_log('Set up done!!!\n')
        if reoptimize:
            now = datetime.datetime.now()
            now = str(now)[:-digit]
            self.write_log(f'[{now}] Reoptmizing current geometry ...\n')
            # To start with the geometries of reactants that are near local minima ... 
            relaxing_path = self.calculator.relax_geometry(molecule,constraints,chg,multiplicity,'Reactant',None,1000) 
            if len(relaxing_path) < 1:
                self.write_log(f'[{datetime.datetime.now()}] Calculation did not terminate properly !!!\n')
                file_directory = os.path.join(self.calculator.error_directory,f'{self.calculator.name}.err')
                self.write_log(f'Check file {file_directory} !!! \n')
                self.write_log(f'Terminating pyMCD ... \n')
                #self.make_restart_file(molecule,constraints,num_steps)
                return []            
            try:
                self.num_force_calls += len(relaxing_path)
            except:
                a = 1
            now = datetime.datetime.now()
            now = str(now)[:-digit]
            self.write_log(f'[{now}] Relaxation finished! Start scanning ....\n')
        else:
            self.write_log('Calculating the energy of current geometry ...\n')
            energy = self.calculator.get_energy(molecule,chg,multiplicity)
            if energy is None:
                self.write_log(f'[{datetime.datetime.now()}] Calculation did not terminate properly !!!\n')
                file_directory = os.path.join(self.calculator.error_directory,f'{self.calculator.name}.err')
                self.write_log(f'Check file {file_directory} !!! \n')
                self.write_log(f'Terminating pyMCD ... \n')
                return []
            molecule.energy = energy
        copied_molecule = molecule.copy()
        copied_molecule.energy = molecule.energy
        trajectory = [copied_molecule] # Trajectory is list of coordinates
        #self.print_status(molecule,constraints)
        energy_list = [molecule.energy]
        if restart:
            self.update_pathway(copied_molecule,'a')
        else:
            self.update_pathway(copied_molecule,'w')
        for iteration in range(num_scan):
            starttime = datetime.datetime.now()
            #self.write_log(f'[{starttime}] Performing one step search ... ({len(trajectory)-1}/{total_num_scans}) completed ... \n')
            # Update/Relax geometry
            self.print_status(molecule,constraints)
            self.update_geometry(molecule,constraints,updates,chg,multiplicity)
            tmp_st = str(datetime.datetime.now())[:-digit]
            content = f'[{tmp_st}] Progress ({iteration+1}/{num_scan}): '
            constraint_contents = []
            for coordinate in constraints:
                new_constraint = self.get_corrected_indices(coordinate)
                constraint_contents.append(f'{new_constraint}: {round(molecule.get_internal_coordinate(coordinate,"degree"),5)}')
            content = content + ', '.join(constraint_contents)
            self.write_log(content+'\n')
            relaxing_path = self.calculator.relax_geometry(molecule,list_of_constraints,chg,multiplicity,'test',self.num_relaxation,self.step_size)
            #scfenergies = calculated_data.scfenergies[]
            #print ('relaxation energy: ',scfenergies[0] - molecule.energy)
            if len(relaxing_path) < 1:
                self.write_log(f'[{datetime.datetime.now()}] Calculation did not terminate properly !!!\n')
                self.write_log(f'Check file {self.calculator.name}.err ...\n')
                self.write_log(f'Terminating pyMCD ... \n')
                self.make_restart_file(molecule,constraints,num_steps)
                return []

            try:
                self.num_force_calls += len(relaxing_path)
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
            endtime = str(endtime)[:-digit]
            delta_time = str(delta_time)[:-digit]
            self.write_log(f'[{endtime}] {x}E{exponent} {self.energy_unit} has {word}. {delta_time} Taken ... Total {self.num_force_calls} force calls performed !!!\n')
            energy_list.append(energy)
            # Save/Copy new molecule
            self.update_pathway(copied_molecule,'a')
            trajectory.append(copied_molecule)
            delta_energy = energy_list[-1] - energy_list[-2]
            print ('Energy increase: ',delta_energy)
            print ('Energy: ',energy_list)
        self.write_log(f'[{datetime.datetime.now()}] Scan completed ...\n')
        self.write_log(f'Total {self.num_force_calls} force calculations performed ...\n')
        self.write_log(f'Total {self.num_hessian_calls} hessian calculations performed ...\n')
        self.write_profile(trajectory)
        et = datetime.datetime.now()
        self.write_log(f'End time: {et}\n')
        self.write_log(f'Taken time: {et-st}\n')
        self.calculator.error_directory = None
        self.step_size = original_step_size 
        self.previous_gradient = None
        self.previous_hessian = None
        self.previous_geometry = None

        return trajectory



