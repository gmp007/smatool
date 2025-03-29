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
from process_stress import process_md_stresses, process_dft_stresses_qe
from calculate_stress import convert_stress_units
from read_write import read_and_write_kpoints,read_options_from_input,update_qe_object
from modify_incar import ChangeDir,WildCard
from shear_tensile_constraint import ShearTensileConstraint
from calculate_stress import apply_indentation_strain



options = read_options_from_input()
dim = options.get("dimensional", "3D")

strain_direction, slip_direction = options.get("slip_parameters").split(", ")
indenter_radius = options.get("indent_radius", 2.0) 

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

# Voigt notation: [epsilon_xx, epsilon_yy, epsilon_zz, gamma_yz, gamma_xz, gamma_xy]
    mask_mapping = {
        'Tensile_x': [1, 0, 0, 0, 0, 0],
        'Tensile_y': [0, 1, 0, 0, 0, 0],
        'Tensile_z': [0, 0, 1, 0, 0, 0],
        'Shear_xy':  [0, 0, 0, 0, 0, 1],
        'Shear_yz':  [0, 0, 0, 1, 0, 0],
        'Shear':  [0, 0, 0, 0, 1, 0], #Shear_xz
        'Tensile_biaxial': [1, 1, 0, 0, 0, 0], #Tensile_biaxial_xy
        'Hydrostatic': [1, 1, 1, 0, 0, 0],
        'indent_strength': [0, 0, 1, 0, 0, 0],  # 
        'ideal_strength': [1, 1, 1, 0, 0, 0],   # 
        'xx_yy_xy': [1, 1, 0, 0, 0, 1]
    }

    
    # Adjust masks based on dimensionality
    if dim == "2D":
        for key in mask_mapping:
            mask_mapping[key][2] = 0  # Zero out epsilon_zz
            mask_mapping[key][3] = 0  # Zero out gamma_yz
            mask_mapping[key][4] = 0  # Zero out gamma_xz

    if dim == "1D":
        for key in mask_mapping:
            mask_mapping[key][1] = 0  # Zero out epsilon_yy
            mask_mapping[key][2] = 0  # Zero out epsilon_zz
            mask_mapping[key][3] = 0  # Zero out gamma_yz
            mask_mapping[key][4] = 0  # Zero out gamma_xz
            mask_mapping[key][5] = 0  # Zero out gamma_xy
            
            
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
           
        if stress_component_list == 'indent_strength':
            # Handle 'indent_strength' separately
            original, cell = apply_indentation_strain(atoms, indenter_radius=indenter_radius, epsilon_0=strain, k=2, dimensionality=dim)
            strained_atoms.set_cell(cell, scale_atoms=True)
        else:
            # Retrieve the mask for the specified stress component
            mask = mask_mapping[stress_component_list]

            # Initialize strain components in Voigt notation
            # [epsilon_xx, epsilon_yy, epsilon_zz, gamma_yz, gamma_xz, gamma_xy]
            strain_components = [strain * m for m in mask]

            # Extract individual strain components
            epsilon_xx = strain_components[0]
            epsilon_yy = strain_components[1]
            epsilon_zz = strain_components[2]
            gamma_yz   = strain_components[3]
            gamma_xz   = strain_components[4]
            gamma_xy   = strain_components[5]

            # Convert engineering shear strains to tensor shear strains
            epsilon_yz = gamma_yz / 2
            epsilon_xz = gamma_xz / 2
            epsilon_xy = gamma_xy / 2

            # Construct the strain tensor
            epsilon_tensor = np.array([
                [epsilon_xx, epsilon_xy, epsilon_xz],
                [epsilon_xy, epsilon_yy, epsilon_yz],
                [epsilon_xz, epsilon_yz, epsilon_zz]
            ])

            # Apply the strain tensor to the cell
            new_cell = np.dot((np.eye(3) + epsilon_tensor), cell)
            strained_atoms.set_cell(new_cell, scale_atoms=True)
        
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
                        energy_convergence_threshold = 1e-4
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



