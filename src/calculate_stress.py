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
import logging
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from ase import Atoms, units
from ase.io import read, write
from ase.cell import Cell
from ase.calculators.vasp import Vasp
from ase.constraints import StrainFilter
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.langevin import Langevin
from ase.io.trajectory import Trajectory
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from read_write import read_options_from_input, write_incar, read_incars, read_and_write_kpoints,load_structure,convert_cif_to_vasp,modify_incar_and_restart
from process_stress import process_md_stresses, process_dft_stresses
from optimize_struct import string_to_tuple

from modify_incar import IncarModifier, ChangeDir,WildCard
from shear_tensile_constraint import ShearTensileConstraint


options = read_options_from_input()
dim = options.get("dimensional", "3D")

strain_direction, slip_direction = options.get("slip_parameters").split(", ")
indenter_radius = options.get("indent_radius", 2.0) 

slip_direction = string_to_tuple(slip_direction, dim)
strain_direction = string_to_tuple(strain_direction, dim)

def convert_stress_units(stress_dict, atoms, dim, valid_components):
    """
    
    Parameters:
    - stress_dict: A dictionary containing average stress components in kbar.
    - atoms: Atoms object with optimized positions.
    - valid_components: A list of valid stress components to process.
    - dim: Dimension of the material ("2D" or "3D").
    
    Returns:
    - A dictionary containing processed stress components in N/m.
    """

    def stress_value(key):
        if dim == "2D":
            length_angl = atoms.cell.cellpar()
            lz_angstrom = length_angl[2]
    
            stress_kbar = stress_dict[key]
            stress_GPa = stress_kbar / 10  # 1 kbar = 0.1 GPa
            converted_stress_value = -1.0 * stress_GPa * lz_angstrom / 10.0
        else:
            # Assuming similar conversion for 3D, modify if needed
            converted_stress_value = -1.0 * stress_dict[key] * 0.10
        return converted_stress_value

    # Creating the dictionary for the stress values using valid components
    converted_stress = {}
    
    for comp in valid_components:
        converted_stress[comp] = stress_value(comp)

    return converted_stress

   
def extract_stress_from_vasprun(xml_path):
    # Parse the XML file
    tree = ET.parse(xml_path)
    root = tree.getroot()
    
    stress_matrix_list = []
    
    # Iterate over each 'varray' with the name 'stress'
    for varray in root.findall(".//varray[@name='stress']"):
        matrix = []
        for v in varray.findall('v'):
            # Split the string into a list of floats
            row = list(map(float, v.text.split()))
            matrix.append(row)
        stress_matrix_list.append(matrix)
    
    return stress_matrix_list



def apply_indentation_strain(atoms, indenter_radius, epsilon_0, k=2.0, dimensionality="2D"):
    """
    Applies radial strain to an atomic structure simulating indentation.

    Args:
        atoms (Atoms): ASE Atoms object representing the structure.
        indenter_radius (float): Radius of the indenter.
        epsilon_0 (float): Maximum strain value at the center of indentation.
        k (float, optional): Decay constant for strain away from the center. Defaults to 2.0.
        dimensionality (str, optional): Dimensionality of the system ("1D", "2D", "3D"). Defaults to "2D".
    """
    
    def radial_strain(r, indenter_radius, epsilon_0, k):
        """Calculates radial strain based on distance from indentation center."""
        return np.where(r < indenter_radius,
                        epsilon_0 * (1 - r / indenter_radius),
                        epsilon_0 * np.exp(-k * (r - indenter_radius)))

    positions = atoms.get_positions()
    center = np.mean(positions, axis=0)
    distances = np.linalg.norm(positions - center, axis=1)
    strains = radial_strain(distances, indenter_radius, epsilon_0, k)

    # Effective Strain Calculation
    strain_threshold = 0.01  # 1%
    indices = np.where(strains < epsilon_0 * strain_threshold)[0] #np.argmax(strains < epsilon_0 * strain_threshold)
    
    if indices.size > 0:
        effective_radius_index = indices[0]
        effective_radius = distances[effective_radius_index]
    else:
        # Set to maximum distance if no strains are below the threshold
        effective_radius = np.max(distances)
        
    #effective_radius = distances[effective_radius_index]

    # Calculate maximum effective strain based on indenter radius
    max_effective_strain_factor = 1.0 + (indenter_radius - effective_radius) / indenter_radius
    effective_strain = epsilon_0 * max_effective_strain_factor
    
    #print("Max Effective Strain Factor:", max_effective_strain_factor)
    #print("Effective Strain:", effective_strain)

    # Update atomic positions based on the calculated strain
    new_positions = positions.copy()
    for i, strain in enumerate(strains):
        displacement = (positions[i] - center) * strain
        new_positions[i] = positions[i] + displacement
    
    atoms.set_positions(new_positions)
    
    # Update cell dimensions based on the calculated strain
    cell = atoms.get_cell()
    new_cell = cell.copy()

    # Plot strain distribution
    plt.figure()
    plt.plot(distances, strains, 'o')
    plt.xlabel('Distance from Center (Å)')
    plt.ylabel('Strain')
    plt.title('Strain Distribution under Indentation')
    plt.grid(True)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    file_to_save = f'strain_distribution_for_{epsilon_0}.png'
    plt.savefig(file_to_save, dpi=400, bbox_inches='tight')
    plt.close()  

    # Construct strain tensor for indentation
    if dimensionality in ["3D", "1D"]:
        epsilon_tensor = np.array([
            [0,          0,          0],
            [0,          0,          0],
            [0,          0, effective_strain]
        ])
    elif dimensionality == "2D":
        epsilon_tensor = np.array([
            [0,          0,          0],
            [0, effective_strain,    0],
            [0,          0,          0]
        ])

    # Apply the strain tensor to the cell
    new_cell = np.dot((np.eye(3) + epsilon_tensor), cell)

#    if dimensionality in ["3D", "1D"]:
#        # Apply normal strain in the z-direction (assumes indentation along the z-axis)
#        new_cell[2, 2] += cell[2, 2] * effective_strain  # Normal strain along the z-axis

#        # Apply shear components if required by the model
#        new_cell[0, 2] += new_cell[0, 0] * effective_strain  # Shear in the z-direction based on the x-component of the first lattice vector
#        new_cell[1, 2] += new_cell[1, 0] * effective_strain  # Shear in the z-direction based on the x-component of the second lattice vector
#        new_cell[2, 2] += new_cell[2, 0] * effective_strain  # Shear in the z-direction based on the x-component of the third lattice vector
    
#    elif dimensionality == "2D":            
#        # Apply normal strain in the y-direction (assuming indentation is in the y-direction)
#        new_cell[1, 1] += cell[1, 1] * effective_strain  # Normal strain along the y-axis

#        # Apply shear components if required by the model
#        new_cell[0, 1] += new_cell[0, 0] * effective_strain  # Shear in the y-direction based on the x-component of the first lattice vector
#        new_cell[1, 0] += new_cell[1, 1] * effective_strain  # Shear in the x-direction based on the y-component of the second lattice vector
  
    atoms.set_cell(new_cell, scale_atoms=False)  # Update the cell of atoms
    
    return cell, new_cell




def pad_to_three_dimensions(vector):
    """
    Pads a 1D or 2D vector with zeros to make it a 3D vector.

    Args:
        vector (np.ndarray): The input vector (1D or 2D).

    Returns:
        np.ndarray: The 3D padded vector.
    """
    return np.pad(vector, (0, 3 - len(vector)), 'constant')

    
    


#slip_direction = pad_to_three_dimensions(np.array(slip_direction))
#strain_direction = pad_to_three_dimensions(np.array(strain_direction))
  

def calculate_stress_for_strain(atoms, strain, stress_component_list=None, dim="2D", mode="DFT"):
    """
    Calculate stress for a given strain using the provided atomic structure.
    
    Parameters:
    - atoms: The ASE Atoms object (the atomic structure).
    - strain: The applied strain.
    - stress_component: Stress component to return ('xx', 'yy', 'xy', etc.). If None, it reads from smatool.in.
    - mode: The mode of calculation, either "DFT" or "MD".

    Returns:
    - The stress for the given strain for the designated component.
    """

#[xx, yy, zz, yz, xz, xy]
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
    strained_atoms = copy.deepcopy(atoms) #read("OPT/optimized_structure.traj").copy()
    cwd = os.getcwd()
    if mode == "DFT":
        kpts, _ = read_and_write_kpoints('static', fileName="KPOINTS-sd", outputDirectory='Yield')
        parent_dir = os.path.dirname(cwd)
        write_incar('yield', cwd, output_dir='Yield')
        incar_settings = read_incars("yield", "INCAR", "Yield")
    elif mode == "MD":
        kpts, _ = read_and_write_kpoints('dynamic', fileName="KPOINTS-sd", outputDirectory='Yield')
        parent_dir = os.path.dirname(cwd)        
        write_incar('mdyield', cwd, output_dir='Yield')
        incar_settings = read_incars("mdyield", "INCAR", "Yield")

    modifyincar = IncarModifier()
    modify_stress = WildCard(incar_settings)
 
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
                temperature_damping_timescale = 100 * units.fs
            
                dt = incar_settings["potim"]    
                total_timesteps = int(options.get("md_timestep", 20))
                modify_stress.modify('nsw', 0)
                temperature_K = incar_settings["tebeg"]
                strained_atoms.set_calculator(Vasp(xc='PBE', kpts=kpts, **incar_settings))
                #dyn = VelocityVerlet(strained_atoms, timestep=dt*units.fs)
                sf = StrainFilter(strained_atoms, mask=mask)
                dyn = NVTBerendsen(strained_atoms, timestep=dt*units.fs, temperature_K=temperature_K,taut=temperature_damping_timescale,fixcm=True)
                dyn.attach(Trajectory('strained_MD.traj', 'w', sf.atoms).write, interval=10)
                MaxwellBoltzmannDistribution(strained_atoms, temperature_K=temperature_K,force_temp=True)

                for _ in range(total_timesteps):
                    constraint = ShearTensileConstraint(stress_component_list, dim)
                    strained_atoms.set_constraint(constraint)
                    dyn.run(steps=1)  
                    strained_atoms.set_cell(sf.get_cell(), scale_atoms=True)
                    #print(f"Step {_}: Stress Temp = {strained_atoms.get_temperature() :.2f} K")
                print("Performing additional thermal stress-strain calculation")
                write('POSCAR', dyn.atoms, format='vasp', direct=True)
                modify_stress.reset()
                os.system(options.get("job_submit_command"))
                #run_monitor_md()
                stress_components = extract_stress_from_vasprun('vasprun.xml')

                converted_stress,_ = process_md_stresses(stress_components, strained_atoms, dim,stress_component_list,portion_to_average=0.5)
             
            elif mode == "DFT":

                try:
                    previous_energy = None
                    energy_converged = False
                    modify_stress.modify('nsw',0)
                    strained_atoms.set_calculator(Vasp(xc='PBE', kpts=kpts, **incar_settings))
                    iteration_counter = 0
                    energy_convergence_threshold = 1e-8  # Starting threshold

                    while not energy_converged:
                        iteration_counter += 1
                        constraint = ShearTensileConstraint(stress_component_list, dim)
                        strained_atoms.set_constraint(constraint)
                        sf = StrainFilter(strained_atoms, mask=mask) #ExpCellFilter(strained_atoms, mask=mask)
                        sf.get_stress()
                        current_energy = sf.get_potential_energy()
                        if previous_energy is not None and abs(current_energy - previous_energy) < energy_convergence_threshold:
                            energy_converged = True

                        # Update the convergence threshold after a certain number of iterations
                        if iteration_counter == 5:
                            energy_convergence_threshold = 1e-4
                        elif iteration_counter == 20:
                            print("Convergence not achieved after 20 iterations. Exiting loop.")
                            break
                        previous_energy = current_energy
                    modify_stress.reset()
                    #shutil.copy('CONTCAR','POSCAR')
                    write("POSCAR", sf.atoms, format='vasp', direct=True)
                    os.system(options.get("job_submit_command"))
                    stress_components = extract_stress_from_vasprun('vasprun.xml')
                    
                    converted_stress = process_dft_stresses(stress_components, strained_atoms, dim, stress_component_list)
                except Exception as e_ase:
                    print("Error occurred in DFT stress-strain:", e_ase)
                    modifyincar.modify_incar(params={'EDIFFG': '-0.09'})  #modify_incar(reset=False) # Check calc.set(prec='Accurate',ediff=1e-2)
                    os.system(options.get("job_submit_command"))
                    stress_components = extract_stress_from_vasprun('vasprun.xml')
                    
                    converted_stress = process_dft_stresses(stress_components, strained_atoms, dim, stress_component_list)
                modifyincar.modify_incar(reset=True)
        except BaseException as e:
            problematicstrain = f"{strain}_vasprun.xml"
            os.system(f"cp vasprun.xml {problematicstrain}") 
            print(f"Error during primary stress calculations: {e}")
            print("Attempting to manually extract stress from vasprun.xml...")
            stress_matrix = extract_stress_from_vasprun(problematicstrain)

            if mode == "MD":
                stresses_xx = [matrix[0][0] for matrix in stress_matrix]
                stresses_yy = [matrix[1][1] for matrix in stress_matrix]
                stresses_zz = [matrix[2][2] for matrix in stress_matrix]
                if dim =="2D":
                    stresses_xy = [(matrix[0][1] + matrix[1][0]) / 2 for matrix in stress_matrix]  
                elif dim =="3D" or dim =="1D":
                    stresses_xy = [(matrix[0][2] + matrix[2][0]) / 2 for matrix in stress_matrix]  

                if dim =="2D":
                    stresses_xz = [(matrix[0][1] + matrix[1][0]) / 2 for matrix in stress_matrix] 
                elif dim =="3D" or dim =="1D":
                    stresses_xz = [(matrix[0][2] + matrix[2][0]) / 2 for matrix in stress_matrix]                   
                
                stresses_xx_yy = [(matrix[0][0] + matrix[1][1]) / 2 for matrix in stress_matrix]
                stresses_xx_yy_xy = [(matrix[0][0] + matrix[1][1] + matrix[0][1] + matrix[1][0] ) / 4 for matrix in stress_matrix]
               
                portion_to_average = 0.5 
                start_idx = int(len(stresses_xx) * (1 - portion_to_average))
                stress_arrays = {
                    'Tensile_x': np.mean(stresses_xx[start_idx:]),
                    'Tensile_y': np.mean(stresses_yy[start_idx:]),
                    'Tensile_z': np.mean(stresses_zz[start_idx:]),
                    'Shear': np.mean(stresses_xy[start_idx:]),
                    'indent_strength': np.mean(stresses_xy[start_idx:]),  
                    'Tensile_biaxial': np.mean(stresses_xx_yy[start_idx:]),
                    'xx_yy_xy': np.mean(stresses_xx_yy_xy[start_idx:])
                }
                
                converted_stress = convert_stress_units(stress_arrays, atoms,dim, stress_component_list)
            elif mode == "DFT":
                stresses_xx = [matrix[0][0] for matrix in stress_matrix]
                stresses_zz = [matrix[2][2] for matrix in stress_matrix]
                stresses_yy = [matrix[1][1] for matrix in stress_matrix]
                if dim == "2D":
                    stresses_xy = [(matrix[0][1] + matrix[1][0]) / 2 for matrix in stress_matrix]
                elif dim =="3D" or dim =="1D":
                    stresses_xy = [(matrix[0][2] + matrix[2][0]) / 2 for matrix in stress_matrix]  

                if dim == "2D":
                    stresses_xz = [(matrix[0][1] + matrix[1][0]) / 2 for matrix in stress_matrix]
                elif dim =="3D" or dim =="1D":
                    stresses_xz = [(matrix[0][2] + matrix[2][0]) / 2 for matrix in stress_matrix]              
                
                stresses_xx_yy = [(matrix[0][0] + matrix[1][1]) / 2 for matrix in stress_matrix]
                stresses_xx_yy_xy = [(matrix[0][0] + matrix[1][1] + matrix[0][1] + matrix[1][0] ) / 4 for matrix in stress_matrix]
                
                stress_arrays = {
                    'Tensile_x': stresses_xx,
                    'Tensile_y': stresses_yy,
                    'Tensile_z': stresses_zz,
                    'Shear': stresses_xy,
                    'indent_strength': stresses_xz,  
                    'Tensile_biaxial': stresses_xx_yy,
                    'xx_yy_xy': stresses_xx_yy_xy
                }
                                
                converted_stress = convert_stress_units(stress_arrays, atoms, dim,stress_component_list)


    if isinstance(converted_stress, (list, np.ndarray)):
        return {component: converted_stress[i] for i, component in enumerate(stress_component_list)}
    else:
        # Handle the scalar case or raise an informative error
        return converted_stress
