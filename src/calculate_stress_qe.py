"""
  SMATool -- Automated toolkit for computing zero and finite-temperature strength of materials

  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation
  version 3 of the License.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE.  See the GNU General Public License for more details.
  Email: cekuma1@gmail.com

""" 

import os
import shutil
import copy
import numpy as np
from ase.calculators.espresso import Espresso
from ase.optimize import LBFGS
from ase.io import read, write
from ase import Atoms, units
from ase.constraints import StrainFilter
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.langevin import Langevin
from ase.io.trajectory import Trajectory
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
#from read_write import applystrain
from optimize_struct_qe import find_qe_pseudopotentials
from process_stress import compute_average_stress, process_md_stresses, process_dft_stresses_qe
from calculate_stress import convert_stress_units
from read_write import read_and_write_kpoints,read_options_from_input,update_qe_object
from modify_incar import ChangeDir,WildCard
from shear_tensile_constraint import ShearTensileConstraint



options = read_options_from_input()

def extract_stress_from_qe_output(output_path):
    stress_matrix_list = []

    # Flag to indicate if the next lines contain stress information
    read_stress = False

    with open(output_path, 'r') as file:
        for line in file:
            if 'total   stress' in line:
                read_stress = True
                continue
            if read_stress:
                if line.strip():  # Non-empty line
                    row = list(map(float, line.split()[3:6]))  
                    stress_matrix_list.append(row)
                else:
                    read_stress = False

    return stress_matrix_list



def sum_of_last_temperatures(temperatures, last_n_runs):
    """
    Calculates the sum of the temperatures of the last n runs.

    Parameters:
    temperatures (list): List of temperatures from the MD simulation.
    last_n_runs (int): Number of last runs to consider for the temperature sum.

    Returns:
    float: Sum of the temperatures of the last n runs, or None if not enough data.
    """
    if len(temperatures) >= last_n_runs:
        return sum(temperatures[-last_n_runs:])
    else:
        print("Not enough data to calculate the temperature sum for the specified number of runs.")
        return None



def calculate_stress_for_strain_qe(atoms, strain, stress_component_list=None, dim="2D", mode="DFT"):
    """
    Calculate stress for a given strain using the provided atomic structure with QE.
    
    Parameters:
    - atoms: The ASE Atoms object (the atomic structure).
    - strain: The applied strain.
    - stress_component_list: Stress component to return.
    - dim: Dimension of the system ("2D" or "3D").
    - mode: The mode of calculation, either "DFT" or "MD".

    Returns:
    - The stress for the given strain for the designated component.
    """

    mask_mapping = {
        'Tensile_x': [0, 1, 1, 0, 0, 0],
        'Tensile_y': [1, 0, 1, 0, 0, 0],
        'Tensile_z': [1, 1, 0, 0, 0, 0],
        'Shear': [0, 0, 0, 0, 1, 0],
        'indent_strength': [0, 0, 0, 0, 1, 0],
        'Tensile_biaxial': [0, 0, 1, 0, 0, 0],
        'xx_yy_xy': [1, 1, 0, 0, 0, 1]
    }
    
    if dim == "2D":
        mask_mapping['Shear'] = [0, 0, 0, 0, 0, 1]
        mask_mapping['indent_strength'] = [0, 0, 0, 0, 0, 1]
    
    
    mask = mask_mapping[stress_component_list] 
    
    
    strained_atoms = copy.deepcopy(atoms) 


    custom_options = options['custom_options']
    base_path = custom_options.get('potential_dir', "./")
    os.environ["VASP_PP_PATH"] = os.path.abspath("./potentials")
    pseudopotentials = find_qe_pseudopotentials(atoms, base_path=base_path)
    

    if mode == "DFT":
        kpts, _ = read_and_write_kpoints('static', fileName="KPOINTS-sd", outputDirectory='Yield')
    elif mode == "MD":
        kpts, _ = read_and_write_kpoints('dynamic', fileName="KPOINTS-sd", outputDirectory='Yield')
        
        
    if dim =="2D":
        qe_parameters = {
            'input_data': {
                'control': {
                    'calculation': 'relax',
                    'restart_mode': 'from_scratch',
                    'pseudo_dir': base_path,
                    'tstress': True,
                    'tprnfor': True,
                    'forc_conv_thr': 0.001,
                    'outdir': './OPT'
                },
                'system': {
                    'ecutwfc': 50,
                    'ecutrho': 500,
                    
                },           
                'electrons': {
                    'conv_thr': 1e-8
                },
              'cell': {
                  'cell_dofree': '2Dshape',
                  'press' : 0.0,
                  'press_conv_thr' : 0.5
              },
              'ions': {},  # Initialize 'ions' key here
            },
            'pseudopotentials': pseudopotentials,
            'kpts': kpts,
        }
    else:
        qe_parameters = {
            'input_data': {
                'control': {
                    'calculation': 'relax',
                    'restart_mode': 'from_scratch',
                    'pseudo_dir': base_path,
                    'tstress': True,
                    'tprnfor': True,
                    'forc_conv_thr': 0.001,
                    'outdir': './OPT'
                },
                'system': {
                    'ecutwfc': 50,
                    'ecutrho': 500,
                    
                },
                'electrons': {
                    'conv_thr': 1e-8
                },
              'cell': {
                  'press' : 0.0,
                  'press_conv_thr' : 0.5
              },
              'ions': {},  # Initialize 'ions' key here
            },
            'pseudopotentials': pseudopotentials,
            'kpts': kpts,
        }
            
    with ChangeDir("Yield") as directory:
        cell = strained_atoms.get_cell().copy()

        strain_factor = (1.0 + strain)
        if stress_component_list == 'Tensile_x':        
            cell[0, 0] *= strain_factor  # Apply strain to x-component of the first lattice vector
            cell[1, 0] *= strain_factor  # Apply strain to x-component of the second lattice vector
            if dim=="3D" or dim=="1D":
                cell[2, 0] *= strain_factor  # Apply strain to x-component of the third lattice vector
        elif stress_component_list == 'Tensile_y':
            cell[0, 1] *= strain_factor  # Apply strain to y-component of the first lattice vector
            cell[1, 1] *= strain_factor  # Apply strain to y-component of the second lattice vector
            if dim=="3D" or dim=="1D":
                cell[2, 1] *= strain_factor  # Apply strain to y-component of the third lattice vector       
#            cell[1, 1] *= (1 + strain)
        elif stress_component_list == 'Tensile_z':
            cell[0, 2] *= strain_factor  # Apply strain to z-component of the first lattice vector
            cell[1, 2] *= strain_factor  # Apply strain to z-component of the second lattice vector
            if dim=="3D" or dim=="1D":
                cell[2, 2] *= strain_factor  # Apply strain to z-component of the third lattice vector
        elif stress_component_list == 'Shear':
            if dim =="3D" or dim == "1D":
                #new_values = [[cell[i, 0] + cell[i, 2] * strain, cell[i, 2] + cell[i, 0] * strain] for i in range(3)]
                
                #for i in range(3):
                #    cell[i, 0], cell[i, 2] = new_values[i]
                
                cell[0, 0] += cell[0, 2] * strain  # Shear in x influenced by z
                cell[0, 2] += cell[0, 0] * strain  # Shear in z influenced by x
                cell[1, 0] += cell[1, 2] * strain  # Shear in x influenced by z
                cell[1, 2] += cell[1, 0] * strain  # Shear in z influenced by x
                cell[2, 0] += cell[2, 2] * strain  # Shear in x influenced by z
                cell[2, 2] += cell[2, 0] * strain  # Shear in z influenced by x
            elif dim == "2D":
                #new_values = [[cell[i, 0] + cell[i, 1] * strain, cell[i, 1] + cell[i, 0] * strain] for i in range(2)]
                #for i in range(2):
                #    cell[i, 0], cell[i, 1] = new_values[i]
                cell[0, 0] += cell[0, 1] * strain  # Shear in x influenced by y
                cell[0, 1] += cell[0, 0] * strain  # Shear in y influenced by x
                cell[1, 0] += cell[1, 1] * strain  # Shear in x influenced by y
                cell[1, 1] += cell[1, 0] * strain  # Shear in y influenced by x

        elif stress_component_list == 'indent_strength': #indent_strength
            if dim =="3D" or dim == "1D":
                cell[0, 2] += cell[0, 0] * strain  # Shear in the z-direction based on the x-component of the first lattice vector
                cell[1, 2] += cell[1, 0] * strain  # Shear in the z-direction based on the x-component of the second lattice vector
                cell[2, 2] += cell[2, 0] * strain  # Shear in the z-direction based on the x-component of the third lattice vector
            elif dim == "2D":            
                # Apply shear strain in the xy-plane
                cell[0, 1] += cell[0, 0] * strain  # Shear in y influenced by x
                cell[1, 1] += cell[1, 0] * strain  # Shear in y influenced by x
                            
        elif stress_component_list == 'Tensile_biaxial':
            # Apply biaxial strain in the xx and yy directions
            cell[0, 0] *= strain_factor  # Strain in the x-direction (xx)
            cell[1, 1] *= strain_factor  # Strain in the y-direction (yy)
            cell[1, 0] *= strain_factor  # Strain in the x-direction (xx) for the second lattice vector
            cell[0, 1] *= strain_factor  # Strain in the y-direction (yy) for the first lattice vector
            if dim=="3D" or dim=="1D":
                cell[2, 0] *= strain_factor  # Strain in the x-direction (xx) for the third lattice vector
                cell[2, 1] *= strain_factor  # Strain in the y-direction (yy) for the third lattice vector
        elif stress_component_list == 'xx_yy_xy':
            # Biaxial strain
            cell[0, 0] *= (1 + strain)  # Strain in the x-direction (xx)
            cell[1, 1] *= (1 + strain)  # Strain in the y-direction (yy)
            # Shear strain
            cell[0, 1] += cell[1, 1] * strain_factor  # Shear in the y-direction based on the y-component of the second lattice vector
            cell[1, 0] += cell[0, 0] * strain_factor  # Shear in the x-direction based on the x-component of the first lattice vector


            # Additionally, apply the strain to the x and y components of all lattice vectors
            cell[1, 0] *= (1 + strain)  # Strain in the x-direction (xx) for the second lattice vector
            cell[0, 1] *= (1 + strain)  # Strain in the y-direction (yy) for the first lattice vector
            if dim=="3D" or dim=="1D":
                cell[2, 0] *= (1 + strain)  # Strain in the x-direction (xx) for the third lattice vector
                cell[2, 1] *= (1 + strain)  # Strain in the y-direction (yy) for the third lattice vector

        strained_atoms.set_cell(cell, scale_atoms=True)
        
        try:
            if mode == "MD":
                md_parameters = qe_parameters.copy()
              #  md_parameters['input_data']['ions'].update({  
              #      'tempw': temperature_K,
              #      'tolp': 1.0e-4
              #  })

                if 'ions' not in qe_parameters['input_data']:
                    qe_parameters['input_data']['ions'] = {}
                
                md_parameters['input_data']['cell'].update({
                    'cell_dynamics': 'bfgs',
                    'press': 0.0,
                    'wmass': 0.02
                })
                md_parameters['input_data']['control']['calculation'] = 'md'   

                            
                md_parameters = update_qe_object("MD Yield Strength", md_parameters)
                md_parameters['input_data']['control'].update({
                    'nstep': 1,
                }) # Ensure nstep 1
                
                ions_params = md_parameters.get('ions', {})
                temperature_K = ions_params.get('temperature',300)
                dt = ions_params.get('dt', 41.341)  # Default time step
                dt = dt/41.341 #Atomic unit to fs for ASE
            
                total_timesteps = int(options.get("md_timestep", 20))
                

                gamma = 0.01/ units.fs  # Friction coefficient
                temperature_damping_timescale = 100 * units.fs
                        
 
                strained_atoms.set_calculator(Espresso(**md_parameters))
                

                sf = StrainFilter(strained_atoms, mask=mask)
                dyn = NVTBerendsen(strained_atoms, timestep=dt*units.fs, temperature_K=temperature_K,taut=temperature_damping_timescale,fixcm=True)
                dyn.attach(Trajectory('strained_MD.traj', 'w', sf.atoms).write, interval=10)
                MaxwellBoltzmannDistribution(strained_atoms, temperature_K=temperature_K,force_temp=True)

                stress_matrices_initial = extract_stress_from_qe_output('../OPT/espresso.pwo')
                c = [stress_matrices_initial]
                temperatures = []
                for _ in range(total_timesteps):
                    constraint = ShearTensileConstraint(stress_component_list, dim)
                    strained_atoms.set_constraint(constraint)
                    dyn.run(steps=1)  
                    strained_atoms.set_cell(sf.get_cell(), scale_atoms=True)
                    #print(f"Step {_}: Stress Temp = {strained_atoms.get_temperature():.2f} K")
                    current_stress_matrix = extract_stress_from_qe_output('espresso.pwo')
                    c.append(current_stress_matrix)
                    temperatures.append(strained_atoms.get_temperature())
                temperature_sum = sum_of_last_temperatures(temperatures, 3)
                #print(temperature_sum)

                stress_components = c
                converted_stress,_ = process_md_stresses(stress_components, strained_atoms, dim,stress_component_list,portion_to_average=0.5)

            elif mode == "DFT":
                temperature_sum = None # Placeholder to enable return tuple
                qe_parameters = update_qe_object("DFT Yield Strength", qe_parameters)
                strained_atoms.set_calculator(Espresso(**qe_parameters))

                previous_energy = None
                energy_converged = False
                iteration_counter = 0
                energy_convergence_threshold = 1e-8  # Starting threshold

                while not energy_converged:
                    iteration_counter += 1
                    constraint = ShearTensileConstraint(stress_component_list, dim)
                    strained_atoms.set_constraint(constraint)
                    sf = StrainFilter(strained_atoms, mask=mask) 
                    sf.get_stress()
                    current_energy = sf.get_potential_energy()
                    if previous_energy is not None and abs(current_energy - previous_energy) < energy_convergence_threshold:
                        energy_converged = True

                    if iteration_counter == 5:
                        energy_convergence_threshold = 1e-6
                    elif iteration_counter == 20:
                        print("Convergence not achieved after 20 iterations. Exiting loop.")
                        break
                    previous_energy = current_energy
                stress_components = extract_stress_from_qe_output('espresso.pwo')
                converted_stress = process_dft_stresses_qe(stress_components, strained_atoms, dim, stress_component_list)

                

        except Exception as e:
            print("Error occurred in QE stress-strain calculation:", e)

    if isinstance(converted_stress, (list, np.ndarray)):
        return {component: converted_stress[i] for i, component in enumerate(stress_component_list)},temperature_sum
    else:
        return converted_stress,temperature_sum



