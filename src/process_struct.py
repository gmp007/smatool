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
 
import sys
import os
from ase.io import read, write

def save_contcar_data(filename, strain, stress, contcar_data):
    if not hasattr(save_contcar_data, "_has_been_called"):
        save_contcar_data._has_been_called = True
        
        if os.path.exists(filename):
            os.remove(filename)

    with open(filename, 'a') as file:
        file.write(f"Strain: {strain}, Stress: {stress}\n")
        file.write(contcar_data + "\n\n")


def get_contcar_data(filename='./Yield/CONTCAR'):
    with open(filename, 'r') as file:
        return file.read()
        
        
        
        
        
def get_contcar_data_from_qe(qe_output_file='./Yield/espresso.pwo'):
    # Read the structure from QE output file
    atoms = read(qe_output_file, format='espresso-out')
    return atoms

def save_contcar_data_qe(filename, strain, stress, atoms):
    if not hasattr(save_contcar_data, "_has_been_called"):
        save_contcar_data._has_been_called = True
        if os.path.exists(filename):
            os.remove(filename)

    with open(filename, 'a') as file:
        file.write(f"Strain: {strain}, Stress: {stress}\n")
        write(file, atoms, format='vasp', vasp5=True, direct=True)
        file.write("\n\n")
        
        
        

