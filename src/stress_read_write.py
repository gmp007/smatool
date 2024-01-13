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
import re



def save_stress_data(filename, component, method, strains, stresses, volumes, energies, pressures, temps=None):
    with open(filename, 'a') as f:
        f.write(f"Component: {component}\n")

        if method == "DFT":
            f.write("Strain, Stress, Mean Volume, Mean Energy, Mean External Pressure \n")
            f.write("0.0, 0.0, 0.0, 0.0 0.0 \n")
            for strain, stress, volume, energy, pressure in zip(strains, stresses, volumes, energies, pressures):
                f.write(f"{strain}, {stress}, {volume}, {energy}, {pressure}\n")

        elif method == "MD":
            f.write("Strain, Stress, Mean Volume, Mean Energy, Mean External Pressure, Mean Temp \n")
            f.write("0.0, 0.0, 0.0, 0.0, 0.0 \n")
            for strain, stress, volume, energy, pressure, temp in zip(strains, stresses, volumes, energies, pressures, temps or []):
                f.write(f"{strain}, {stress}, {volume}, {energy},  {pressure}, {temp}\n")
                
        f.write("\n")
        
        
        
def read_stress_data_2nd(filename):
    """
    Reads stress data from a file and returns the data as a dictionary of lists.

    Parameters:
    filename (str): The name of the file containing the data.

    Returns:
    dict: A dictionary with keys 'strains', 'stresses', 'volumes', 'energies', 'temps', and 'pressures'.
    """
    strain, stress, volume, energy, temperature, pressure = [], [], [], [], [], []

    with open(filename, 'r') as file:
        lines = file.readlines()[1:]  # Skip the header

    for line in lines:
        try:
            values = line.split()
            strain.append(float(values[0]))
            stress.append(float(values[1]))  # Stress in GPa
            volume.append(float(values[2]))
            energy.append(float(values[3]))
            temperature.append(float(values[4]))
            pressure.append(float(values[5]))
        except:
            pass  # Skip lines that can't be parsed

    return {
        "strains": strain,
        "stresses": stress,
        "volumes": volume,
        "energies": energy,
        "temps": temperature,
        "pressures": pressure
    }




        
def read_stress_data(filename, component, file_type, method):
    strains, stresses, volumes, energies, pressures, temps = [], [], [], [], [], []
    reading = False

    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(f"Component: {component}"):
                reading = True
                continue

            if reading:
                if line.strip() == "" or "Strain" in line:
                    continue  # Skip empty lines and headers

                data = line.split(',')
                if file_type == "VASP":
                    if len(data) >= 5:  # VASP usually has at least 5 data columns
                        if data[1].strip().lower() == 'none':
                            continue
                        try:
                            strains.append(float(data[0].strip()))
                            stresses.append(float(data[1].strip()))
                            volumes.append(float(data[2].strip()))
                            energies.append(float(data[3].strip()))
                            pressures.append(float(data[4].strip()))
                            if method == "MD" and len(data) >= 6:
                                temps.append(float(data[5].strip()))
                        except ValueError:
                            continue  # Skip lines that cannot be converted to float
                elif file_type == "QE":
                    if len(data) >= 4:  # QE usually has at least 4 data columns
                        try:
                            strains.append(float(data[0].strip()))
                            stresses.append(float(data[1].strip()))
                            volumes.append(float(data[2].strip()))
                            energies.append(float(data[3].strip()))
                            if method == "MD" and len(data) >= 5:
                                temps.append(float(data[4].strip()))
                        except ValueError:
                            continue  # Skip lines that cannot be converted to float

    return {
        "strains": strains,
        "stresses": stresses,
        "volumes": volumes,
        "energies": energies,
        "pressures": pressures, # if file_type == "VASP" else [],
        "temps": temps if method == "MD" else []
    }


# Usage
# result = read_stress_data('path/to/outputfile', 'desired_component', 'VASP', 'MD')


        
        

def mean_value(file_name, num_last_samples, keyword, column_name):
    values = []

    with open(file_name, 'r') as file:
        # Read lines from the file
        lines = file.readlines()

        for line in lines:
            if keyword in line:
                if column_name == 'temperature':
                    #match = re.search(r'mean temperature.*:\s+(\d+\.\d+)', line) #re.search(r'\d+\.\d+', line)
                    match = re.search(r'mean temperature.*:\s+(\d+\.\d+)', line)
                    if match:
                        values.append(float(match.group(1)))

                elif column_name == 'volume':
                    match = re.search(r'\d+\.\d+', line)
                    if match:
                        values.append(float(match.group()))
                elif column_name == 'pressure':
                    match = re.search(r'external pressure =\s+(-?[\d.]+)', line)
                    if match:
                        values.append(float(match.group(1)))
                elif column_name == 'totalenergy':
                    match = re.search(r'free  energy   TOTEN  =\s+(-?[\d.]+)', line)
                    if match:
                        values.append(float(match.group(1)))                       
    if len(values) >= num_last_samples:
        return sum(values[-num_last_samples:]) / num_last_samples
    elif values:
        return sum(values) / len(values)
    else:
        return None

def wildcard_output(outcar_file, vasprun_file, num_last_samples):
    mean_temp = mean_value(outcar_file, num_last_samples, 'mean temperature', 'temperature')
    mean_vol = mean_value(vasprun_file, num_last_samples, '<i name="volume">', 'volume') #You can also use OUTCAR "volume of cell"
    mean_press = mean_value(outcar_file, num_last_samples, 'external pressure =', 'pressure')
    total_energy = mean_value(outcar_file, num_last_samples, 'free  energy   TOTEN  =', 'totalenergy')  

    return mean_temp, mean_vol, total_energy, -mean_press/10.




#QE mainly
def mean_value_qe(file_name, num_last_samples, regex_pattern):
    values = []

    with open(file_name, 'r') as file:
        for line in file:
            match = re.search(regex_pattern, line)
            if match:
                try:
                    value = float(match.group(1))
                    values.append(value)
                except ValueError:
                    pass  # Ignore lines where conversion to float fails

    # Compute the mean of the last 'num_last_samples' values
    if len(values) >= num_last_samples:
        return sum(values[-num_last_samples:]) / num_last_samples
    elif values:
        return sum(values) / len(values)
    else:
        return None

def wildcard_output_qe(pw_out_file, num_last_samples):
    energy_regex = r'!\s+total energy\s+=\s+(-?[\d.]+)\s+Ry'  # For total energy
    volume_regex = r'unit-cell volume\s+=\s+([\d.]+)\s+\(a\.u\.\)\^3'  # For volume in (a.u.)^3
    temperature_regex = r'temperature\s+=\s+([\d.]+)\s+K'  # For temperature
    pressure_regex = r'P=\s+(-?[\d.]+)'  # For hydrostatic pressure

    mean_energy = mean_value_qe(pw_out_file, num_last_samples, energy_regex)
    mean_vol_au = mean_value_qe(pw_out_file, num_last_samples, volume_regex)
    mean_temp = mean_value_qe(pw_out_file, num_last_samples, temperature_regex)
    mean_press = mean_value_qe(pw_out_file, num_last_samples, pressure_regex)

    # Convert mean energy from Rydberg to electronvolts
    if mean_energy is not None:
        mean_energy *= 13.6056980659

    # Convert mean volume from (a.u.)^3 to Angstrom^3
    if mean_vol_au is not None:
        mean_vol_angstrom = mean_vol_au * (0.529177**3)

    return mean_temp, mean_vol_angstrom, mean_energy, -mean_press/10.
    


