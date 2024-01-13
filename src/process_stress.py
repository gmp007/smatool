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

import numpy as np
from get_youngs_modulus import get_youngs_modulus
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

def compute_average_stress(stress_component_list, dim,stress_components, atoms):
    if dim == "2D":
        length_angl = atoms.cell.cellpar()
        lz_angstrom = length_angl[2]
        conversion_lz = -1.0*lz_angstrom / 100.
    else:
        conversion_lz = -0.10
    
    if stress_component_list == "Tensile_x":
        #print("This is the stress ", stress_components[0, 0] * conversion_lz)
        return  stress_components[0, 0] * conversion_lz
    elif stress_component_list == "Tensile_y":
        return  stress_components[1, 1] * conversion_lz
    elif stress_component_list == "Tensile_z":
        return  stress_components[2, 2] * conversion_lz
    elif stress_component_list == "Shear":
        if dim =="2D":
            return stress_components[0, 1] * conversion_lz
        elif dim =="3D" or dim =="1D":
            return stress_components[0, 2] * conversion_lz
    elif stress_component_list == "indent_strength":
        if dim =="2D":
            return stress_components[0, 1] * conversion_lz
        elif dim =="3D" or dim =="1D":
            return stress_components[0, 2] * conversion_lz
    elif stress_component_list == "Tensile_biaxial":
        return  (stress_components[0, 0] + stress_components[1, 1]) * conversion_lz/2.
    elif stress_component_list == "xx_yy_xy":
        return  (stress_components[0, 0] + stress_components[1, 1] + stress_components[0, 1])/3. * conversion_lz
    else:
        raise ValueError(f"Unknown stress component: {stress_component_list}")
        

def process_md_stresses(stresses, atoms, dim,stress_component_list, portion_to_average=0.5):
    if not isinstance(stresses, np.ndarray):
        stresses = np.array(stresses)
    #print("This is the shape ", stresses.shape)
    dim1, _, _ = stresses.shape
    #print(stresses, stress_component_list)



    stress_components = np.zeros((dim1, 3, 3))

    # Aggregate the stress components into shape (5,3,3)
    #print("stresses ", stresses)
    for i in range(dim1):
        stress_matrix = stresses[i]
        stress_components[i, 0, 0] = stress_matrix[0, 0]  # xx component
        stress_components[i, 1, 1] = stress_matrix[1, 1]  # yy component
        stress_components[i, 2, 2] = stress_matrix[2, 2]  # zz component
        stress_components[i, 0, 1] = stress_matrix[0, 1]  # xy component
        stress_components[i, 1, 0] = stress_matrix[1, 0]  # yx component (symmetrical)
        stress_components[i, 0, 2] = stress_matrix[0, 2]  # xz component
        stress_components[i, 2, 0] = stress_matrix[2, 0]  # zx component (symmetrical)
        stress_components[i, 1, 2] = stress_matrix[1, 2]  # yz component
        stress_components[i, 2, 1] = stress_matrix[2, 1]  # zy component (symmetrical)

    # Calculate average stress over runs
    average_stresses = np.zeros((3, 3))
    start_idx = int(dim1 * (1 - portion_to_average))
    for i in range(3):
        for j in range(3):
            average_stresses[i, j] = np.mean(stress_components[start_idx:, i, j])

    # Get the last stress value
    last_stress = stress_components[-1]
    def map_stress_values(stress, stress_component_list):
        mapped_values = []
        for component in stress_component_list.split(','):
            value = compute_average_stress(component, dim,stress, atoms)
            mapped_values.append(value)
        return mapped_values[0]
    
    average_mapped_stress = map_stress_values(average_stresses, stress_component_list)
    last_mapped_stress = map_stress_values(last_stress, stress_component_list)
    return average_mapped_stress, last_mapped_stress



def process_dft_stresses_qe(stresses, atoms, dim, stress_component_list):
    if not isinstance(stresses, np.ndarray):
        stresses = np.array(stresses)

    # Check if stresses array has the correct shape (9x3)
    if stresses.shape[1] != 3 or stresses.shape[0] % 3 != 0:
        raise ValueError("Stresses array does not have the correct shape.")

    # Extract and reshape the last three rows into a 3x3 matrix
    last_stress_matrix = stresses[-3:].reshape((3, 3))

    # Initialize the stress components array
    stress_components = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            stress_components[i, j] = last_stress_matrix[i][j]

    # Symmetrize the off-diagonal components
    for i in range(3):
        for j in range(i + 1, 3):
            symmetrical_stress = (stress_components[i, j] + stress_components[j, i]) / 2
            stress_components[i, j] = stress_components[j, i] = symmetrical_stress

    # Map the stress components to the requested values
    mapped_stress_values = [compute_average_stress(component, dim, stress_components, atoms)
                            for component in stress_component_list.split(',')]

    return mapped_stress_values[0]





def process_dft_stresses(stresses, atoms, dim, stress_component_list):
    # Convert the list to a numpy array if it isn't one already
    if not isinstance(stresses, np.ndarray):
        stresses = np.array(stresses)
    
    # Extract the stress matrix from the last run
    last_stress_matrix = stresses[-1]

    # Process the components
    stress_components = np.zeros((3, 3))
    stress_components[0, 0] = last_stress_matrix[0][0]  # xx component
    stress_components[1, 1] = last_stress_matrix[1][1]  # yy component
    stress_components[2, 2] = last_stress_matrix[2][2]  # zz component
    stress_components[0, 1] = (last_stress_matrix[0][1] + last_stress_matrix[1][0]) / 2  # xy component
    stress_components[0, 2] = (last_stress_matrix[0][2] + last_stress_matrix[2][0]) / 2  # xz component
    stress_components[1, 0] = stress_components[0, 1]  # yx component
    stress_components[1, 2] = (last_stress_matrix[1][2] + last_stress_matrix[2][1]) / 2  # yz component
    stress_components[2, 0] = stress_components[0, 2]  # zx component
    stress_components[2, 1] = stress_components[1, 2]  # zy component

    # Map the stress components
    mapped_stress_values = []
    for component in stress_component_list.split(','):
        value = compute_average_stress(component, dim, stress_components, atoms)
        mapped_stress_values.append(value)
    mapped_stress_values = mapped_stress_values[0]
    
    return mapped_stress_values
    
   

def get_yield_strength_offsetold(stresses, strains, offset=0.002, tolerance=1e-4):
    """
    Calculates the yield strength of a material given stresses and strains.

    :param stresses: List of stress values
    :param strains: List of strain values corresponding to the stresses
    :param offset: Offset for the yield strength calculation (default is 0.002)
    :param tolerance: Numerical tolerance for comparing stress values
    :return: Calculated yield strength
    """

    stresses = np.array(stresses, dtype=float)
    strains = np.array(strains, dtype=float)

    min_length = min(len(stresses), len(strains))
    stresses = stresses[:min_length]
    strains = strains[:min_length]
  
    
    E = get_youngs_modulus(stresses, strains)
    
    if E is None:
        raise ValueError("Young's modulus could not be determined.")

    interpolation = interp1d(strains, stresses, kind='linear', bounds_error=False, fill_value="extrapolate")

    min_strain, max_strain = np.min(strains), np.max(strains)
    interpolated_strains = np.linspace(min_strain, max_strain, num=50)
    interpolated_stresses = interpolation(interpolated_strains)
    
    
    for stress, strain in zip(interpolated_stresses, interpolated_strains):
        offset_stress = E * (strain - offset)
        if stress >= offset_stress - tolerance:  # Allow for some numerical tolerance
            yield_strength = stress
            break

    if yield_strength <= 0:
        print("Warning: Yield strength is non-positive. Returning maximum stress as heuristic yield strength.")
        return np.max(stresses)

    return yield_strength




def get_yield_strength_offset(stresses, strains, offset=0.002, tolerance=1e-4):
    stresses = np.array(stresses, dtype=float)
    strains = np.array(strains, dtype=float)

    min_length = min(len(stresses), len(strains))
    stresses = stresses[:min_length]
    strains = strains[:min_length]
    
    if len(stresses) != len(strains):
        raise ValueError("Stresses and strains arrays must be of the same length.")

    E = get_youngs_modulus(stresses, strains)

    if E is None:
        raise ValueError("Young's modulus could not be determined.")

    # Interpolating the stress-strain curve
    interpolation = interp1d(strains, stresses, kind='linear', bounds_error=False, fill_value="extrapolate")

    # Calculate stress at zero strain (extrapolation)
    stress_at_zero_strain = interpolation(0)
    # Searching for the yield strength
    for strain in np.linspace(np.min(strains), np.max(strains), num=1000):
        offset_stress = E * offset + stress_at_zero_strain
        actual_stress = interpolation(strain)
        if actual_stress >= offset_stress - tolerance:
            #print("actual_stress  ", actual_stress)
            #exit(0)
            return actual_stress

    raise ValueError("Yield strength could not be determined within the given strain range.")
    


def get_yield_strength_linearity(strain, stress):
    """
    Determines yield strength finding deviation from linear regime of a stress-strain curve.

    Parameters:
    strain (numpy array): An array of strain values.
    stress (numpy array): An array of corresponding stress values.

    Returns:
    tuple: A tuple containing the maximum deviation from linearity,
           the strain at which this maximum deviation occurs, and
           the parameters of the best linear fit.
    """

    def linear_fit(x, a, b):
        return a * x + b
    params, _ = curve_fit(linear_fit, strain, stress)

    extended_strain_range = np.linspace(strain.min(), strain.max(), num=100)
    extended_linear_fit_values = linear_fit(extended_strain_range, *params)

    deviations = stress - np.interp(strain, extended_strain_range, extended_linear_fit_values)
    max_deviation = stress[np.argmax(np.abs(deviations))] #np.max(np.abs(deviations))
    max_deviation_strain = strain[np.argmax(np.abs(deviations))]

    return max_deviation, max_deviation_strain, params
