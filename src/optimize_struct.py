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
import numpy as np
import copy
import spglib
from ase import Atoms, units
from ase.io import read, write
from ase.calculators.vasp import Vasp
from ase.geometry import cell_to_cellpar, cellpar_to_cell
from ase.spacegroup import get_spacegroup, crystal
from ase.optimize import LBFGS, BFGS
from ase.constraints import StrainFilter
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.io.trajectory import Trajectory
from shear_tensile_constraint import ShearTensileConstraint
from modify_incar import IncarModifier, ChangeDir,WildCard
from read_write import read_options_from_input,write_incar, read_incars, read_and_write_kpoints,load_structure,modify_incar_and_restart, simplify_formula
from euler_angles_rotation import apply_rotation, rotate, rotate_crystal_structure, apply_rotation_2D, rotate_crystal_structure_2D
from euler_angles_rotation import calculate_schmid_factor, generate_slip_systems_2d_hexagonal, generate_slip_systems


options = read_options_from_input()
slipon = options.get("slipon",False)
checkrotation = options.get("rotation", False)

#checkrotation = False
#if slipon:
#    rotation = True
#    checkrotation = rotation
#else:
#    checkrotation = options.get("rotation", False)
   
#checkrotation = options.get("rotation",False)
dim = options.get("dimensional", "3D")
use_saved_data = options.get('use_saved_data', False)

def remove_spurious_distortion(pos):
    # Normalize and orthogonalize the cell vectors
    cell_params = cell_to_cellpar(pos.get_cell())
    new_cell = cellpar_to_cell(cell_params)
    pos.set_cell(new_cell, scale_atoms=True)

    # Adjust atom positions
    pos.wrap()

    pos.center()

    return pos
    

def swap_axes_to_longest_c(pos):
    """
    Reorders the lattice vectors of the given structure (pos) such that the longest vector is always c-axis.

    Parameters:
    pos (ASE Atoms object): The atomic structure.

    Returns:
    ASE Atoms object: Updated structure with reordered lattice vectors.
    """
    a, b, c = pos.get_cell()
    lengths = [np.linalg.norm(a), np.linalg.norm(b), np.linalg.norm(c)]
    max_index = lengths.index(max(lengths))

    new_cell = [a, b, c]  
    new_cell[2], new_cell[max_index] = new_cell[max_index], new_cell[2]  

    if max_index == 0:
        new_cell[0], new_cell[1] = new_cell[1], new_cell[0]
    pos.set_cell(new_cell)

    return pos

    
def run_calculation(atoms, calculator_settings, fmax=0.02, max_retries=5, retry_count=0):
    try:
        if atoms.get_calculator() is None:
            atoms.set_calculator(Vasp(**calculator_settings))

        opt = LBFGS(atoms)
        opt.run(fmax=fmax)
        #atoms.get_potential_energy()
    except Exception as e:
        if retry_count < max_retries:
            new_fmax = 0.02 + 0.01 * retry_count
            print(f"Caught an exception: {e}")
            print("Modifying INCAR and restarting the calculation.")
            modify_incar_and_restart()

            atoms = read("POSCAR")  # Re-read atoms from POSCAR
            atoms.set_calculator(Vasp(**calculator_settings))  # Reset calculator
            run_calculation(atoms, calculator_settings, new_fmax, max_retries, retry_count + 1)
            
        else:
            print("Maximum number of retries reached. Exiting.")
    print("DFT Optimization Done!")




def string_to_tuple(s, dim="3D"):
    result = []
    i = 0
    while i < len(s):
        if s[i] == '-':
            # Ensure the next character is a digit and combine it with the minus sign
            if i + 1 < len(s) and s[i + 1].isdigit():
                result.append(-int(s[i + 1]))
                i += 2
        elif s[i].isdigit():
            result.append(int(s[i]))
            i += 1
        else:
            # Skip any non-digit, non-minus characters
            i += 1
    
    if dim == "2D" and len(result) > 2:
        # Remove the last element if dim is "2D" and there are more than 2 elements
        result = result[:2]

    return tuple(result)


def string_to_tupleold(s):
    result = []
    i = 0
    while i < len(s):
        if s[i] == '-':
            # Ensure the next character is a digit and combine it with the minus sign
            if i + 1 < len(s) and s[i + 1].isdigit():
                result.append(-int(s[i + 1]))
                i += 2
        elif s[i].isdigit():
            result.append(int(s[i]))
            i += 1
        else:
            # Skip any non-digit, non-minus characters
            i += 1
    return tuple(result)
    

def check_vasp_optimization_completed(atoms, mode, output_dir="OPT"):
    """
    Checks if the VASP optimization has already been completed.

    Parameters:
    atoms (Atoms): ASE Atoms object of the initial structure.
    mode (str): Mode of the optimization.
    output_dir (str): Directory where VASP output files are stored.

    Returns:
    bool: True if optimization is completed, False otherwise.
    """
    optimized = False
    contcar_file = os.path.join(output_dir, "CONTCAR")
    outcar_file = os.path.join(output_dir, "OUTCAR")

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    if os.path.isfile(contcar_file) and os.path.isfile(outcar_file):
        with open(outcar_file, "r") as f:
            content = f.read()
            if "reached required accuracy " in content:
                optimized_atoms = read(contcar_file)
                if set(atoms.get_chemical_symbols()) == set(optimized_atoms.get_chemical_symbols()):
                    print("DFT Optimization already completed. Skipping...")
                    optimized = True
                    write(os.path.join(output_dir, "CONTCAR"), atoms, format='vasp', direct=True)
                else:
                    print("Structures in CONTCAR and initial atoms object do not match. Proceeding with optimization...")

    return optimized
                
def optimize_structure(mode,stress_component_list):
    """
    Optimizes the given atomic structure based on the specified mode.

    Parameters:
    - atoms: The ASE Atoms object (the atomic structure).
    - mode: Optimization mode (either 'DFT' or 'MD').
    """
    from ase import Atoms

        
    struct = options.get("structure_file")
    atoms = load_structure(struct)
    if dim =="3D":
        cell = atoms.get_cell()
        a_length = np.linalg.norm(cell[0])
        b_length = np.linalg.norm(cell[1])
        c_length = np.linalg.norm(cell[2])
        if max(a_length, b_length, c_length) > 12.5: #Gauge for nanoribbon
            atoms = swap_axes_to_longest_c(atoms)
    atoms = remove_spurious_distortion(atoms)
    spg = get_spacegroup(atoms)


#    if "P" in spg.symbol:
#        # The structure is primitive; convert to conventional
#        conventional_cell = crystal(atoms.get_chemical_symbols(), 
#                                    basis=atoms.get_scaled_positions(), 
#                                    spacegroup=spg.no, 
#                                    cellpar=atoms.get_cell_lengths_and_angles())
#        atoms = conventional_cell

#

    #rotation_dir = options.get("strain_direction")
    strain_direction, slip_direction = options.get("slip_parameters").split(", ")
    slip_direction = string_to_tuple(slip_direction, dim)
    strain_direction = string_to_tuple(strain_direction, dim)
        
    if not os.path.exists("OPT"):
        os.mkdir("OPT")
        
    cwd = os.getcwd()
    if mode == "DFT":
        write_incar('opt', cwd, output_dir='OPT')
        kpts, _ = read_and_write_kpoints('static', fileName="KPOINTS-sd", outputDirectory='OPT')
        incar_settings = read_incars("opt", "INCAR", "OPT")
        atoms.set_calculator(Vasp(xc='PBE', kpts=kpts, **incar_settings))

    elif mode == "MD":
        write_incar('opt', cwd, output_dir='OPT')
        kpts, _ = read_and_write_kpoints('static', fileName="KPOINTS-sd", outputDirectory='OPT')
        incar_settings = read_incars("opt", "INCAR", "OPT")
        atoms.set_calculator(Vasp(xc='PBE', kpts=kpts, **incar_settings))

        kpts_nostress, _ = read_and_write_kpoints('dynamic', fileName="KPOINTS-sd", outputDirectory='OPT')
        write_incar('mdnostress', cwd, output_dir='OPT')
        incar_settings_nostress = read_incars("mdnostress", "INCAR", "OPT")

        
    else:
        raise ValueError("Invalid mode specified. Choose either 'DFT' or 'MD'.")

    # Perform structure optimization
    optimized = check_vasp_optimization_completed(atoms, mode)
    with ChangeDir("OPT"):
        atoms.set_calculator(Vasp(xc='PBE', kpts=kpts, **incar_settings))
        calculator_settings = {'xc': 'PBE', 'kpts': kpts, **incar_settings}

        if not optimized and not use_saved_data:        
            run_calculation(atoms,calculator_settings)
            optimized = True

        if checkrotation and optimized:
            print("Atoms object updated with with optimized structure ")
            #rotation_dir = options.get("strain_direction")
            write("current_structure.vasp", atoms,format='vasp', direct=True)
            if dim =="2D":
                rotate_crystal_structure_2D("current_structure.vasp", strain_direction, slip_direction, dim, atoms,slipon)
            else:            
                rotate_crystal_structure("current_structure.vasp", strain_direction, slip_direction, dim, atoms,slipon)
            atoms = read("POSCAR_rotated")
            write("structure_rotated.cif", atoms)
            #os.system("rm -rf POSCAR_rotated")
            os.system("rm -rf current_structure.vasp")
            
            # Calculate Schmid factor
            strain_dir = strain_direction #string_to_tuple(options.get("strain_direction"))
            if dim =="2D":
                slip_systems = generate_slip_systems_2d_hexagonal()
                schmid_factors = [calculate_schmid_factor(normal, direction, strain_dir) for normal, direction in slip_systems]
                try:
                    lattice = atoms.cell.array
                    positions = atoms.get_scaled_positions()
                    numbers = atoms.numbers
                    cell = (lattice, positions, numbers)                   
                    crystal_system = spglib.get_spacegroup(cell,symprec=0.1)
                except TypeError:
                    crystal_system = spglib.get_spacegroup(atoms,symprec=0.1)
            elif dim =="3D":
                slip_systems,crystal_system = generate_slip_systems(atoms)
                schmid_factors = [calculate_schmid_factor(normal, direction, strain_dir) for normal, direction in slip_systems]

            finite_unique_schmid_factors = set(factor for factor in schmid_factors if np.isfinite(factor))
            formatted_schmid_factors = ", ".join("{:.3f}".format(factor) for factor in finite_unique_schmid_factors)
            chem_formula = atoms.get_chemical_formula(mode='hill',empirical=False)
            chem_formula = simplify_formula(chem_formula)
            print(f"Unique set of Schmid factor for {crystal_system} {chem_formula} in {strain_dir} strain direction are: \n ",  '(',formatted_schmid_factors,')' )
            file_path = os.path.join("..", "schmid_factor.dat")
            np.savetxt(file_path, np.array(list(finite_unique_schmid_factors)), fmt='%.3f', header=f'Unique set of Schmid factor for {crystal_system} {chem_formula} in {strain_dir} strain direction are:')
            print("Schmid factors saved to file.")
    

        if mode == "MD":
        #[xx, yy, zz, yz, xz, xy]
            
            mask_mapping = {
                'Tensile_x': [1, 0, 0, 0, 0, 0],
                'Tensile_y': [0, 1, 0, 0, 0, 0],
                'Tensile_z': [0, 0, 1, 0, 0, 0],
                'Shear': [0, 0, 0, 0, 1, 0],
                'indent_strength': [0, 0, 0, 0, 1, 0],
                'Tensile_biaxial': [1, 1, 0, 0, 0, 0],
                'xx_yy_xy': [1, 1, 0, 0, 0, 1]
            }
            
            if dim == "2D":
                mask_mapping['Shear'] = [0, 0, 0, 0, 0, 1]
                mask_mapping['indent_strength'] = [0, 0, 0, 0, 0, 1]
                
        
            mask = mask_mapping[stress_component_list] 
            modify_nostress = WildCard(incar_settings_nostress)

            mdsteps = int(options.get("md_timestep",20))
            modify_nostress.modify('nsw',0)
            supercell_dimensions = options.get('md_supercell', [1, 1, 1])

            temperature_K = incar_settings_nostress["tebeg"]
            gamma = 0.01/units.fs  # Friction coefficient
            dt = incar_settings_nostress["potim"]
            atoms = atoms.repeat(supercell_dimensions)
            atoms.set_calculator(Vasp(xc='PBE', kpts=kpts_nostress, **incar_settings_nostress))
            sf = StrainFilter(atoms,mask=mask)
            pre_dyn = Langevin(atoms, timestep=dt*units.fs, temperature_K=temperature_K, friction=gamma,fixcm=True)
            pre_dyn.attach(Trajectory('equilibration.traj', 'w', sf.atoms).write, interval=10) 
            MaxwellBoltzmannDistribution(atoms, temperature_K=temperature_K,force_temp=True)

            for _ in range(mdsteps):   
                constraint = ShearTensileConstraint(stress_component_list, dim)
                atoms.set_constraint(constraint)
                pre_dyn.run(steps=1)
                atoms.set_cell(sf.get_cell(), scale_atoms=True)
                #print(f"Step {_}: Equilibration Temp = {atoms.get_temperature():.2f} K")
                   
                modify_nostress.reset()
            print("Equilibration Done!")

    write("OPT/optimized_structure.traj", atoms)

    dir_yield = "Yield"
    if os.path.exists(dir_yield):
        shutil.rmtree(dir_yield)
    os.mkdir(dir_yield)
    
    
    optimized_atoms = Atoms(symbols=atoms.get_chemical_symbols(), positions=atoms.get_positions(), cell=atoms.get_cell(), pbc=True)
    structure_file = "OPT/optimized_structure.cif"
    write(structure_file, optimized_atoms)
    return optimized_atoms
