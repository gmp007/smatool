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
from ase.geometry import cell_to_cellpar, cellpar_to_cell
from ase.spacegroup import get_spacegroup, crystal
from ase.calculators.espresso import Espresso
from ase.constraints import StrainFilter
from ase.md.langevin import Langevin
from ase.optimize import LBFGS, BFGS
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.io.trajectory import Trajectory
from ase.units import fs, kB
from optimize_struct import swap_axes_to_longest_c, remove_spurious_distortion,string_to_tuple
from shear_tensile_constraint import ShearTensileConstraint
from read_write import read_and_write_kpoints,read_options_from_input,load_structure,update_qe_object,simplify_formula
from euler_angles_rotation import calculate_euler_angles, apply_rotation, rotate, rotate_crystal_structure, apply_rotation_2D, rotate_crystal_structure_2D
from euler_angles_rotation import calculate_schmid_factor, generate_slip_systems_2d_hexagonal, generate_slip_systems
from modify_incar import ChangeDir,WildCard
from pathlib import Path
import json

options = read_options_from_input()
slipon = options.get("slipon",False)
checkrotation = options.get("rotation",False)
dim = options.get("dimensional", "3D")
use_saved_data = options.get('use_saved_data', False)


def find_qe_pseudopotentials(atoms, base_path="./potentials"):
    # Convert relative path to absolute path
    #base_path = os.path.abspath(base_path)
    base_path = Path(base_path).absolute()

    pseudopotentials = {}
    unique_symbols = set(atoms.get_chemical_symbols())

    for symbol in unique_symbols:
        potential_paths = [
            os.path.join(base_path, symbol + "_pdojo.upf"),
            os.path.join(base_path, symbol + ".UPF"),
            os.path.join(base_path, symbol + "pz-vbc.UPF"),
            os.path.join(base_path, symbol + "_sv.UPF"),
            os.path.join(base_path, symbol + ".upf")
        ]

        pp_file_path = next((path for path in potential_paths if os.path.exists(path)), None)

        if not pp_file_path:
            raise Exception(f"Pseudopotential for {symbol} not found in any of the expected directories!")

        pseudopotentials[symbol] = os.path.basename(pp_file_path)


    return pseudopotentials




def check_optimization_completed(atoms, mode, output_dir="OPT"):
    """
    Checks if the DFT optimization has already been completed.

    Parameters:
    atoms (Atoms): ASE Atoms object of the initial structure.
    mode (str): Mode of the optimization.
    output_dir (str): Directory where output files are stored.

    Returns:
    bool: True if optimization is completed, False otherwise.
    """
    optimized = False
    output_file = os.path.join(output_dir, "espresso.pwo")
    structure_file = os.path.join(output_dir, "optimized_structure.traj")  

    if os.path.exists(output_dir):
        if os.path.isfile(output_file):
            with open(output_file, "r") as f:
                content = f.read()
                if "End of self-consistent calculation" in content: 
                    # Check if optimized structure file exists and matches the initial atoms
                    if os.path.isfile(structure_file):
                        optimized_atoms = read(structure_file)
                        if set(atoms.get_chemical_symbols()) == set(optimized_atoms.get_chemical_symbols()):
                            print("DFT Optimization already completed. Skipping...")
                            optimized = True
                        else:
                            print("Structures in final structure file and initial atoms object do not match. Proceeding with optimization...")
                    else:
                        print("Optimized structure file not found. Proceeding with optimization...")
    else:
        os.mkdir(output_dir)

    return optimized


    

def run_calculation(atoms, qe_parameters, fmax=0.02, max_retries=5, retry_count=0):
    try:
        if atoms.get_calculator() is None:
            atoms.set_calculator(Espresso(**qe_parameters))

        opt = LBFGS(atoms)
        opt.run(fmax=fmax)
        # You can also obtain energy, forces, etc., as needed
        # energy = atoms.get_potential_energy()

    except Exception as e:
        if retry_count < max_retries:
            new_fmax = 0.02 + 0.01 * retry_count
            print(f"Caught an exception: {e}")
            print("Modifying QE parameters and restarting the calculation.")

            atoms = read("optimized_structure.traj")  # Re-read atoms from QE output
            atoms.set_calculator(Espresso(**qe_parameters))  # Reset calculator
            run_calculation(atoms, qe_parameters, new_fmax, max_retries, retry_count + 1)
            
        else:
            print("Maximum number of retries reached. Exiting.")
    print("DFT Optimization Done!")
    

def optimize_structure_qe(mode, stress_component_list):
    """
    Optimizes the given atomic structure based on the specified mode using Quantum ESPRESSO.

    Parameters:
    - mode: Optimization mode ('DFT' or 'MD').
    - stress_component_list: List of stress components.
    - options: Dictionary of options including structure file path and other parameters.
    """
    
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


    strain_direction, slip_direction = options.get("slip_parameters").split(", ")
    slip_direction = string_to_tuple(slip_direction)
    strain_direction = string_to_tuple(strain_direction)
    


    custom_options = options.get('custom_options', {})
    base_path = custom_options.get('potential_dir', "./potentials")
    os.environ["VASP_PP_PATH"] = os.path.abspath(base_path)

        
    #base_path = options.get('potential_dir', "./potentials")
    pseudopotentials = find_qe_pseudopotentials(atoms, base_path=base_path)

    optimized = check_optimization_completed(atoms, mode)
    kpts, _ = read_and_write_kpoints('static', fileName="KPOINTS-sd", outputDirectory='OPT')
    #kpts = [2, 2, 1]
    
    if dim =="2D":
        qe_parameters = {
            'input_data': {
                'control': {
                    'calculation': 'vc-relax',
                    'restart_mode': 'from_scratch',
                    'pseudo_dir': base_path,
                    'tstress': True,
                    'tprnfor': True,
                    'forc_conv_thr': 0.001,
                    'outdir': './OPT'
                },
                'system': {
                    'ecutwfc': 50,
                    'ecutrho': 600,
                    'occupations': 'smearing',
                    'smearing': 'mp',
                    'degauss' : 0.02,

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
                    'calculation': 'vc-relax',
                    'restart_mode': 'from_scratch',
                    'pseudo_dir': base_path,
                    'tstress': True,
                    'tprnfor': True,
                    'forc_conv_thr': 0.001,
                    'outdir': './OPT'
                },
                'system': {
                    'ecutwfc': 50,
                    'ecutrho': 600,
                    
                },
                'electrons': {
                    'conv_thr': 1e-8
                },
              'cell': {
                  'cell_dofree': 'all',
                  'press' : 0.0,
                  'press_conv_thr' : 0.5
              },
            'ions': {},  # Initialize 'ions' key here
            },
            'pseudopotentials': pseudopotentials,
            'kpts': kpts,
        }

    with ChangeDir("OPT"):
        qe_parameters = update_qe_object("DFT Optimization", qe_parameters)
        atoms.set_calculator(Espresso(**qe_parameters))
        calculator_settings = qe_parameters
        optimized_atoms = atoms.copy()
        
                
        if not optimized and not use_saved_data:        
            run_calculation(optimized_atoms,calculator_settings)
            optimized = True

        if checkrotation and optimized:
            print("Atoms object updated with optimized structure ")
            write("current_structure.vasp", atoms,format='vasp', direct=True)
            if dim =="2D":
                rotate_crystal_structure_2D("current_structure.vasp", strain_direction, slip_direction, dim, atoms,slipon)
            else:            
                rotate_crystal_structure("current_structure.vasp", strain_direction, slip_direction, dim, atoms,slipon)
            atoms = read("POSCAR_rotated")
            write("structure_rotated.cif", atoms)
            os.system("rm -rf POSCAR_rotated")
            os.system("rm -rf current_structure.vasp")

            # Calculate Schmid factor
            strain_dir = strain_direction #string_to_tuple(options.get("strain_direction"))
            if dim =="2D":
                slip_systems = generate_slip_systems_2d_hexagonal()
                schmid_factors = [calculate_schmid_factor(normal, direction, strain_dir) for normal, direction in slip_systems]
                crystal_system = spglib.get_spacegroup(atoms,symprec=0.1)
            elif dim =="3D":
                slip_systems,crystal_system = generate_slip_systems(atoms)
                schmid_factors = [calculate_schmid_factor(normal, direction, strain_dir) for normal, direction in slip_systems]
            #
            finite_unique_schmid_factors = set(factor for factor in schmid_factors if np.isfinite(factor))
            formatted_schmid_factors = ", ".join("{:.3f}".format(factor) for factor in finite_unique_schmid_factors)
            chem_formula = atoms.get_chemical_formula(mode='hill',empirical=False)
            chem_formula = simplify_formula(chem_formula)
            print(f"Unique set of Schmid factor for {crystal_system} {chem_formula} in {strain_dir} strain direction are: \n ",  '(',formatted_schmid_factors,')' )
            file_path = os.path.join("..", "schmid_factor.dat")
            np.savetxt(file_path, np.array(list(finite_unique_schmid_factors)), fmt='%.3f', header=f'Unique set of Schmid factor for {crystal_system} {chem_formula} in {strain_dir} strain direction are:')
            print("Schmid factors saved to file.")
                
    
        if mode == "MD":    

            mdsteps = int(options.get("md_timestep",20))

            gamma = 0.01/fs  # Friction coefficient
            
            md_parameters = qe_parameters.copy()

            if 'ions' not in qe_parameters['input_data']:
                qe_parameters['input_data']['ions'] = {}

            md_parameters['input_data']['cell'].update({
                'press': 0.0,
                'wmass': 0.02
            })

            
            md_parameters['input_data']['control']['calculation'] = 'vc-md'

            md_parameters['input_data']['control'].update({
                'nstep': 1,
            })
            md_parameters = update_qe_object("MD Nostrain", md_parameters)

            ions_params = md_parameters.get('ions', {})
            temperature_K = ions_params.get('temperature',300)
            dt = ions_params.get('dt', 41.341)  # Default time step
            dt = dt/41.341 #Atomic unit to fs for ASE
            
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


            supercell_dimensions = options.get('md_supercell', [1, 1, 1])
            atoms = atoms.repeat(supercell_dimensions)
                        
            atoms.set_calculator(Espresso(**md_parameters))
            sf = StrainFilter(atoms,mask=mask)
            
            pre_dyn = Langevin(atoms, timestep=dt*fs, temperature_K=temperature_K, friction=gamma,fixcm=True)
            pre_dyn.attach(Trajectory('equilibration.traj', 'w', sf.atoms).write, interval=10) 
            MaxwellBoltzmannDistribution(atoms, temperature_K=temperature_K,force_temp=True)
            

            for _ in range(mdsteps):   
                constraint = ShearTensileConstraint(stress_component_list, dim)
                atoms.set_constraint(constraint)
                pre_dyn.run(steps=1)
                atoms.set_cell(sf.get_cell(), scale_atoms=True)
                temperature = atoms.get_temperature() 
                #print(f"Step {_}: Equilibration Temp = {temperature:.2f} K")
                #print("Forces ", atoms.get_forces())
            print("Equilibration Done!")

    write("OPT/optimized_structure.traj", atoms)

    # Write the optimized structure to a CIF file
    structure_file = "OPT/optimized_structure.cif"
    write(structure_file, optimized_atoms)
        
    dir_yield = "Yield"
    if os.path.exists(dir_yield):
        shutil.rmtree(dir_yield)
    os.mkdir(dir_yield)

    optimized_atoms = Atoms(symbols=atoms.get_chemical_symbols(), positions=atoms.get_positions(), cell=atoms.get_cell(), pbc=True)
    return optimized_atoms


