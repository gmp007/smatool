"""
  THICK-2D -- Thickness Hierarchy Inference & Calculation Kit for 2D materials

  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation
  version 3 of the License.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE.  See the GNU General Public License for more details.
  Email: cekuma1@gmail.com

""" 
 
import sys
import os
import json
from ase.io import read, write
import shutil
import glob
from math import gcd
import re
from write_inputs import write_default_input, write_default_ystool_in, print_default_input_message_0, print_default_input_message_1, print_default_input_message_0
from smatoolgui import display_help, smatoolguicall

def read_options_from_input():

    cwd = os.getcwd()

    gui_flag = False
    if len(sys.argv) > 1:
        first_arg = sys.argv[1].lower()
        gui_flag = first_arg in ["-gui"]


    if gui_flag:
        smatoolguicall()
        sys.exit(0)  
        
        
    ystool_in_exists = os.path.exists(os.path.join(cwd, "smatool.in"))
    run_mode_flag = (len(sys.argv) > 1 and sys.argv[1] == "-0")
    if run_mode_flag and not ystool_in_exists:
        write_default_ystool_in(cwd)
        print_default_input_message_0()
        sys.exit(0)

        
    """
    Read the stress component options from the 'smatool.in' file.
    ...
    rotation = on/off
    ...
    Returns:
    - Dictionary of the component settings.
    """
    options = {
        'code_type': 'VASP',
        'use_saved_data': False,
        'use_dnn_gan': False,
        'strains': [0.0, 0.6, 0.05],
        'yieldpoint_offset': 0.002,
        'nlayers': 1,
        'vdwgap': 3.5,
        'dimensional': [], 
        'mode': 'DFT',
        "scrystal": 'cubic',
        'plusminus': False,
        'use_schmid_factors': False,
        'components': [],
        'md_supercell': [1, 1, 1],
        'md_timestep': 1000,
        'indent_radius': 0.0,
        'custom_options': {},  # to store user-defined options
        'job_submit_command': None,
        'structure_file': None,
        'rotation': False,
        'save_struct': False,
        'slipon': False,
        'slip_parameters': None,
    }

    try:
        with open("smatool.in", "r") as f:
            lines = f.readlines()
           # previous_line = None
            for line in lines:
                line = line.strip()
                if line.startswith("#") or not line:
                    #previous_line = line.strip()
                    continue

                try:
                    key, value = line.split('=', 1)
                    key = key.strip()
                    value = value.strip()
                except ValueError:
                    print(f"Skipping invalid line: '{line}'. Expected format 'key=value'.")
                    continue

                if key in ["structure_file", "job_submit_command", "strains", "slip_parameters"]:
                    options[key] = value
                elif key == "components":
                    options[key] = value.split()
                elif key == "plusminus":
                    val_lower = value.lower()
                    if val_lower in ["on", "true", "yes", "1"]:
                        # Tension + Compression
                        options["plusminus"] = True
                    elif val_lower in ["off", "false", "no", "0"]:
                        # Tension only
                        options["plusminus"] = False
                    elif val_lower == "compress":
                        # Compression only
                        options["plusminus"] = "compress"
                    else:
                        print(f"Warning: Unrecognized plusminus value '{value}'. Defaulting to False.")
                        options["plusminus"] = False
                elif key == "md_supercell":
                    options[key] = [int(i) for i in value.split(',')]
                elif key in ["mode","dimensional","scrystal","code_type"]:
                    options[key] = value.upper()
                elif key in ["use_saved_data","save_struct","rotation", "slipon","use_dnn_gan"]:
                    options[key] = value.lower() in ['true', 'yes', '1','on']
                #elif key in ["rotation", "slipon"]:
                #    options[key] = value.lower() == 'on'
                elif key in options:
                    if key in ["md_timestep", "yieldpoint_offset",'nlayers','vdwgap','indent_radius']:
                        options[key] = float(value)
                    else:
                        options[key] = value.lower() == 'on'
                else:
                    options['custom_options'][key] = value

        # Set the environment variable for ASE's VASP command
        if options.get('job_submit_command'):
            os.environ["ASE_VASP_COMMAND"] = options['job_submit_command']

    except FileNotFoundError:
        print("'smatool.in' file not found. Using default settings.")
        
    dim = options.get("dimensional")
    code_type = options.get("code_type")
    run_mode_flag = (len(sys.argv) > 1 and sys.argv[1] == "-0") #and 'dimensional' in options
    if run_mode_flag and ystool_in_exists:
        write_default_input(cwd, dim,code_type)
        print_default_input_message_1()
        sys.exit(0)


    help_arg = (len(sys.argv) > 1 and (sys.argv[1] == "-help" or sys.argv[1] == "--help" or sys.argv[1] == "--h" or sys.argv[1] == "-h"))
    if help_arg:
        display_help()
        sys.exit(0)

    if options.get("slipon", False):
        options["rotation"] = True
    else:
        options["rotation"] = False
         

    # Validation for code_type and job_submit_command using regular expressions
    job_submit_command = options.get("job_submit_command", "").lower()
    if not options.get("use_saved_data"): 
        if code_type == "VASP" and not re.search(r'vasp(_std)?|vasp_\w+', job_submit_command):
            print("Warning: For code_type 'VASP', please set job_submit_command to the VASP executable (e.g., vasp_std, vasp_*).")
            return None 
        elif code_type == "QE" and not re.search(r'pw\.x', job_submit_command):
            print("Warning: For code_type 'QE', please set job_submit_command to the Quantum Espresso executable (e.g., pw.x).")
            return None
     
               
    return options


    

def write_incar(step, cwd, output_dir=None):
    infile = cwd + '/INCARs'
    outfile = 'INCAR' if output_dir is None else os.path.join(output_dir, 'INCAR')
    tag = '# Step:'
    
    step_dict = {
        'opt': 'DFT Optimization',
        'yield': 'DFT Yield Strength',
        'mdopt': 'MD Optimization',
        'mdnostress': 'MD Nostrain',
        'mdyield': 'MD Yield Strength'
    }

    is_write = False

    #print(f"Looking for: {step_dict[step]}")  # Debug print
    with open(infile, 'r') as infile, open(outfile, 'w') as ofile:
        for line in infile:
            if tag in line:
                # print(f"Found step line: {line.strip()}")  # Debug print
                is_write = step_dict[step] in line
                # print(f"is_write set to: {is_write}")  # Debug print

            if is_write:
                # print(f"Writing line: {line.strip()}")  # Debug print
                ofile.write(line)
   



def read_incars(step, fileName="INCARs", directory=None):
    """
    Reads settings for the specified step and section from the INCARs file.

    Parameters:
    - step: The step corresponding to a section in the INCARs.
    - fileName: The name of the INCARs file. Defaults to "INCARs".
    - directory: The directory where the INCARs file is located.

    Returns:
    - Dictionary containing the settings for the specified step and section.
    """
    step_dict = {
    'opt': 'DFT Optimization',
    'yield': 'DFT Yield Strength',
    'mdopt': 'MD Optimization',
    'mdnostress': 'MD Nostrain',
    'mdyield': 'MD Yield Strength'
    }  
    
    section_name = step_dict[step] 
    
    if directory:
        fileName = os.path.join(directory, fileName)
    
    with open(fileName, "r") as f:
        lines = f.readlines()

    in_section = False
    settings = {}

    for line in lines:
        line = line.strip()

        if line.startswith("#"):
            found_section = line[2:].strip()
            
            if section_name in found_section:
                in_section = True
                continue

        if in_section:
            if "=" in line:
                
                key, value = [item.strip() for item in line.split("=")]
                #print(f"Found step line: {key}, {value}")  # Debug print
                value = value.strip().lower()
                if value in [".true.", "true"]:
                    value = True
                elif value in [".false.", "false"]:
                    value = False
                else:
                    try:
                        value = int(value)
                    except ValueError:
                        try:
                            value = float(value)
                        except ValueError:
                            pass
                settings[key.lower()] = value

    return settings



def read_and_write_kpoints(step=None, calculation_type=None, fileName="KPOINTS", directory=None, outputDirectory=None):
    step_dict = {
        'static': 'Step: Static Calculation',
        'dynamic': 'Step: Dynamical Calculation'
    }
    
    section_name = step_dict.get(step)
    
    if directory:
        fileName = os.path.join(directory, fileName)
    
    with open(fileName, 'r') as f:
        lines = f.readlines()
    
    in_section = False
    kpoints_data = []
    
    for line in lines:
        line = line.strip()
        
        if line.startswith("#"):
            found_section = line[2:].strip()
            in_section = found_section == section_name
        elif in_section:
            kpoints_data.append(line)
    
    # If a specific output directory is given, create the KPOINTS file there
    if outputDirectory:
        output_file_path = os.path.join(outputDirectory, "KPOINTS")
        if not os.path.exists(outputDirectory):
            os.makedirs(outputDirectory)
    else:
        output_file_path = "KPOINTS"
    
    with open(output_file_path, 'w') as f:
        f.write("\n".join(kpoints_data))
    
    kpoints = list(map(int, kpoints_data[2].split()))
    
    return kpoints, output_file_path




def load_structure(options):
    options = read_options_from_input()
    filename = options.get('structure_file', None)
    
    # If a filename is provided in the options, use it directly
    if filename:
        if not os.path.exists(filename):
            raise FileNotFoundError(f"Provided file {filename} is not found.")
        return read(filename)

    # If no filename is provided, proceed with the original logic
    vasp_files = glob.glob("*.vasp")
    cif_files = glob.glob("*.cif")

    if len(vasp_files) > 1 or len(cif_files) > 1:
        raise RuntimeError("Multiple VASP or CIF files detected. Please have only one of each.")

    if vasp_files and cif_files:
        vasp_base = os.path.splitext(vasp_files[0])[0]
        cif_base = os.path.splitext(cif_files[0])[0]

        if vasp_base == cif_base:
            raise RuntimeError("Both VASP and CIF files are present with the same basename. Ambiguous input.")

    if vasp_files:
        return read(vasp_files[0])
    elif cif_files:
        cif_filename = cif_files[0]
        vasp_equivalent = os.path.splitext(cif_filename)[0] + ".vasp"
        cif_orig = os.path.splitext(cif_filename)[0] + "_orig.cif"
        
        convert_cif_to_vasp(cif_filename, vasp_equivalent)  
        os.rename(cif_filename, cif_orig)
        
        return read(vasp_equivalent)
    else:
        raise FileNotFoundError("Neither VASP nor CIF file is found.")



def convert_cif_to_vasp(cif_file, vasp_file):
    """
    Convert a CIF file to VASP format using ASE.
    
    Args:
    - cif_file (str): Path to the input CIF file.
    - vasp_file (str): Path to the output VASP file.
    """
    
    # Read the structure from the CIF file
    atoms = read(cif_file)

    # Write the structure in VASP format
    write(vasp_file, atoms, format='vasp', direct=True)



def modify_incar_and_restart():
    # Decrease EDIFF value in INCAR file
    with open("INCAR", "r") as file:
        lines = file.readlines()

    with open("INCAR", "w") as file:
        for line in lines:
            if line.startswith("EDIFF"):
                # Assuming the line is something like "EDIFF = 1E-4"
                parts = line.split("=")
                ediff_value = float(parts[1].strip())
                new_ediff = ediff_value / 10  # Decrease EDIFF by an order of magnitude
                line = f"EDIFF = {new_ediff}\n"
            file.write(line)

    # Copy CONTCAR to POSCAR for restart
    if os.path.exists("CONTCAR"):
        shutil.copyfile("CONTCAR", "POSCAR")
        


def merge_qe_parameters(existing_params, new_step_data):
    new_input_data = new_step_data.get('input_data', {})
    new_pseudopotentials = new_step_data.get('pseudopotentials', None)
    new_kpts = new_step_data.get('kpts', None)

    original_pseudo_dir = existing_params.get('input_data', {}).get('control', {}).get('pseudo_dir', None)

    for key in ['control', 'system', 'electrons', 'cell', 'ions']:
        if key in new_input_data:
            if key in existing_params.get('input_data', {}):
                if key == 'control':
                    temp_control = new_input_data['control'].copy()
                    temp_control.pop('pseudo_dir', None) 
                    existing_params['input_data']['control'].update(temp_control)
                else:
                    existing_params['input_data'][key].update(new_input_data[key])
            else:
                existing_params['input_data'][key] = new_input_data[key]

    if original_pseudo_dir is not None:
        existing_params['input_data']['control']['pseudo_dir'] = original_pseudo_dir

    if new_pseudopotentials is not None:
        existing_params['pseudopotentials'] = new_pseudopotentials
    if new_kpts is not None:
        existing_params['kpts'] = new_kpts

    return existing_params


def update_qe_object(step_name, existing_qe_parameters, file_name="qe_input.in"):
    if os.path.exists(file_name):
        file_path = file_name
    elif os.path.exists(os.path.join('..', file_name)):
        file_path = os.path.join('..', file_name)
    else:
        raise FileNotFoundError(f"File '{file_name}' not found in the current or parent directory.")

    with open(file_path, 'r') as file:
        data = json.load(file)['steps']

    for step in data:
        if step['name'] == step_name:
            updated_qe_parameters = merge_qe_parameters(existing_qe_parameters, step['data'])
            return updated_qe_parameters

    raise ValueError(f"Step '{step_name}' not found in the file.")


def validate_options(options, valid_options):
    errors = []
    for key, valid_values in valid_options.items():
        value = options.get(key)
        if value not in valid_values:
            errors.append(f"Invalid value for '{key}': {value}. Valid options are: {', '.join(valid_values)}.")
    return errors


def validate_active_components(active_components, valid_components):
    errors = []
    for component in active_components:
        if component not in valid_components:
            errors.append(f"Invalid active component: {component}. Valid options are: {', '.join(valid_components)}.")
    return errors
    
def simplify_formula(formula):
    # Parse the formula into elements and their counts
    elements = re.findall('([A-Z][a-z]*)(\d*)', formula)
    elements = {element: int(count) if count else 1 for element, count in elements}

    # Find the greatest common divisor of the counts
    common_divisor = gcd(*elements.values())

    # Simplify the formula
    simplified_formula = ''.join(f"{element}{(count // common_divisor) if count > common_divisor else ''}" 
                                 for element, count in elements.items())
    
    return simplified_formula
            
if __name__ == '__main__':
    import os
    cwd = os.getcwd()
