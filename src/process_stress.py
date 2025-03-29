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
import pandas as pd
import logging
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.integrate import trapz
from scipy.stats import linregress
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
#from ase.eos import EquationOfState
from ase.units import kJ
    
        
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
        logging.error("Stresses array does not have the correct shape.")

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
    
   



def calculate_yield_strength(strains, stresses, dim, thickness, yieldpoint_offset, filename, start_percent=5, end_percent=95, increment_percent=5):
    """
    Calculates the yield strength from stress-strain data.

    Args:
    strains (array-like): Array of strain values.
    stresses (array-like): Array of stress values.
    dim (str): Dimension of the stress ('2D' or '3D').
    yieldpoint_offset (float): Offset for yield point calculation.
    filename (str): Filename for saving the plot.
    start_percent (int): Starting percentage for segment analysis.
    end_percent (int): Ending percentage for segment analysis.
    increment_percent (int): Increment percentage for segment analysis.

    Returns:
    float: Elastic modulus (E)
    float: Yield strength
    """
    stresses = np.array(stresses, dtype=float)
    strains = np.array(strains, dtype=float)
    
    # Function to find the best linear fit
    def find_best_linear_fit(strains, stresses, segment_percentages):
        best_mse = float('inf')
        best_segment = 0
        for percent in segment_percentages:
            segment = int(len(strains) * percent / 100)
            if segment > 5:  # Ensuring enough points for regression
                segment_strains = strains[:segment]
                segment_stresses = stresses[:segment]
                slope, intercept, r_value, p_value, std_err = linregress(segment_strains, segment_stresses)
                mse = np.mean((segment_stresses - (slope * segment_strains + intercept))**2)
                if mse < best_mse:
                    best_mse = mse
                    best_segment = segment
        return best_segment

    try:

    # Identifying the linear regime
        segment_percentages = range(start_percent, end_percent, increment_percent)
        optimal_segment = find_best_linear_fit(strains, stresses, segment_percentages)

        linear_strains = strains[:optimal_segment] if optimal_segment > 0 else strains
        linear_stresses = stresses[:optimal_segment] if optimal_segment > 0 else stresses

        #linear_strains = strains[:optimal_segment]
        #linear_stresses = stresses[:optimal_segment]
    
        # Calculating the elastic modulus (E) in the identified linear regime
        if len(linear_strains) == 0 or len(linear_stresses) == 0:
            raise ValueError("Linear strains or stresses array is empty.")
            logging.error("Linear strains or stresses array is empty.")
        
        E_al, intercept, _, _, _ = linregress(linear_strains, linear_stresses)

        # Offset Method
        stress_offset_al = E_al * (strains - yieldpoint_offset)

        # Finding the intersection point
        intersection_index = np.argmin(np.abs(stresses - stress_offset_al))
        yield_strength = stresses[intersection_index]
        elastic_strain_at_yield = strains[intersection_index]
        
        # Calculate resilience
        resilience = 0.5 * yield_strength * elastic_strain_at_yield

    except ValueError as e:
        print(f"An error occurred: {e}")
        logging.error(f"An error occurred: {e}")
        E_al = None  
        yield_strength = None  
        resilience = None  
        elastic_strain_at_yield = None

    #Calculate toughness
    toughness = trapz(stresses, strains)
    if dim == "2D":
        resilience = resilience/thickness #*1E-09  
        toughness  = toughness/thickness #*1E-09 

    #Work Hardening Rate
    dσ = np.gradient(np.array(stresses))
    dε = np.gradient(np.array(strains))
    work_hardening_rate = dσ / dε
    work_hardening_rate_df = pd.DataFrame({'#Strain': strains, 'Work Hardening Rate': work_hardening_rate})
    
    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(strains, stresses, label='Stress-Strain Curve')
    plt.plot(strains, stress_offset_al, label='0.2% Offset Line', color='red')
    plt.scatter(strains[intersection_index], yield_strength, color='green', label='Yield Strength')
    plt.xlabel('Strain', fontsize=20)
    ys_percentage = yieldpoint_offset * 100
    ylabel = 'Stress (N/m)' if dim == "2D" else 'Stress (GPa)'
    plt.ylabel(ylabel, fontsize=20)
    plt.title(f'Stress-Strain Curve with {ys_percentage}% Offset Line')
    plt.legend()
    plt.grid(True)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    file_to_save = f'{filename}.png'
    plt.savefig(file_to_save, dpi=400, bbox_inches='tight')
    plt.show()

    return E_al, yield_strength,resilience, toughness, work_hardening_rate_df,elastic_strain_at_yield




def calculate_bulk_modulus(volumes, energies, pressures, filename,dim="2D"):
    """
    Calculates the bulk modulus and its pressure derivative from volume-energy data using the Birch-Murnaghan equation.

    Args:
    volumes (array-like): Array of volume values.
    energies (array-like): Array of energy values.
    pressures (array-like): Array of pressure values.
    filename (str): Filename for saving the plot.

    Returns:
    float: Bulk modulus (K)
    float: Bulk modulus pressure derivative (K')
    """

    def birchmurnaghan(V, E0, B0, BP, V0):
        """
        Birch-Murnaghan equation from PRB 70, 224107.
        """
        eta = (V0 / V)**(1 / 3)
        E = E0 + 9 * B0 * V0 / 16 * (eta**2 - 1)**2 * (
            6 + BP * (eta**2 - 1) - 4 * eta**2)
        return E

    volumes = np.array(volumes, dtype=float)
    energies = np.array(energies, dtype=float)
    pressures = np.array(pressures, dtype=float)

    # Initial guess for the fitting parameters
    p0 = [min(energies), 1, 1, min(volumes)]
    
    try:
        popt, pcov = curve_fit(birchmurnaghan, volumes, energies, p0, maxfev=5000)
        e0, B, Bp, v0 = popt
        
        # Convert bulk modulus from eV/Å^3 to GPa
        B_GPa = B / kJ * 1.0e24
        
        # Plotting the data and the fitted curve
        plt.figure()
        plt.scatter(volumes, energies, label='Data')
        V_fit = np.linspace(min(volumes), max(volumes), 100)
        E_fit = birchmurnaghan(V_fit, *popt)
        plt.plot(V_fit, E_fit, label='Birch-Murnaghan Fit', color='red')

        if dim =="2D":
            plt.text(0.05, 0.95, f'E0: {e0:.3f} eV\nV0: {v0:.3f} Å^2\nB: {B_GPa:.3f} GPa',
                    transform=plt.gca().transAxes,
                    verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))
            plt.xlabel('Volume (Å^3)')
        else:
            plt.text(0.05, 0.95, f'E0: {e0:.3f} eV\nV0: {v0:.3f} Å^3\nB: {B_GPa:.3f} GPa',
                    transform=plt.gca().transAxes,
                    verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))
            plt.xlabel('Volume (Å^3)')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.title('Equation of State Fit')
        plt.savefig(f'{filename}.png')
        plt.show()
        
    except RuntimeError as e:
        print(f"An error occurred: {e}")
        B_GPa, Bp = None, None    

    return B_GPa, Bp
    
    
    
def calculate_bulk_modulus2old(volumes, pressures, filename, dim, start_percent=5, end_percent=95, increment_percent=5):
    """
    Calculates the bulk modulus from pressure-volume data.

    Args:
    volumes (array-like): Array of volume values.
    pressures (array-like): Array of pressure values.
    filename (str): Filename for saving the plot.
    start_percent (int): Starting percentage for segment analysis.
    end_percent (int): Ending percentage for segment analysis.
    increment_percent (int): Increment percentage for segment analysis.

    Returns:
    float: Bulk modulus (K)
    float: Bulk modulus pressure derivative (K')
    """
    volumes = np.array(volumes, dtype=float)
    pressures = np.array(pressures, dtype=float)

    if np.all(volumes == volumes[0]):
        print("Error: All volume values are the same. Cannot calculate bulk modulus.")
        return None, None

    # Function to find the best linear fit
    def find_best_linear_fit(volumes, pressures, segment_percentages):
        best_mse = float('inf')
        best_segment = 0
        for percent in segment_percentages:
            segment = int(len(volumes) * percent / 100)
            if segment > 5:  # Ensuring enough points for regression
                segment_volumes = volumes[:segment]
                segment_pressures = pressures[:segment]
                slope, intercept, r_value, p_value, std_err = linregress(segment_volumes, segment_pressures)
                mse = np.mean((segment_pressures - (slope * segment_volumes + intercept))**2)
                if mse < best_mse:
                    best_mse = mse
                    best_segment = segment
        return best_segment

    try:
        # Identifying the linear regime
        segment_percentages = range(start_percent, end_percent, increment_percent)
        optimal_segment = find_best_linear_fit(volumes, pressures, segment_percentages)

        linear_volumes = volumes[:optimal_segment] if optimal_segment > 0 else volumes
        linear_pressures = pressures[:optimal_segment] if optimal_segment > 0 else pressures

        # Calculating the bulk modulus in the identified linear regime
        if len(linear_volumes) == 0 or len(linear_pressures) == 0:
            raise ValueError("Linear volumes or pressures array is empty.")
        
        slope, intercept, _, _, _ = linregress(linear_volumes, linear_pressures)

        # Calculate the average volume in the elastic region
        average_volume = np.mean(linear_volumes)

        # Calculate the bulk modulus
        bulk_modulus = average_volume * slope

        sorted_indices = np.argsort(volumes)
        volumes_sorted = volumes[sorted_indices]
        pressures_sorted = pressures[sorted_indices]

        # Calculate the derivative of the bulk modulus with respect to pressure (K')
        dp_dv = np.gradient(pressures_sorted, volumes_sorted)
        d2p_dv2 = np.gradient(dp_dv, volumes_sorted)

        # Calculate K' using the second derivative
        bulk_modulus_derivative = average_volume * np.mean(d2p_dv2[:optimal_segment])

    except ValueError as e:
        print(f"An error occurred: {e}")
        bulk_modulus = None
        bulk_modulus_derivative = None

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(volumes, pressures, 'o-', label='Pressure-Volume Curve')
    if optimal_segment > 0:
        plt.plot(linear_volumes, intercept + slope * linear_volumes, 'r--', label='Linear Fit')

    if dim == "2D":
        plt.xlabel('Area (Å²)', fontsize=20)
        plt.ylabel('Pressure (N/m)', fontsize=20)
        plt.title('Pressure-Area Curve with Linear Fit')
    else:
        plt.xlabel('Volume (Å³)', fontsize=20)
        plt.ylabel('Pressure (GPa)', fontsize=20)
        plt.title('Pressure-Volume Curve with Linear Fit')
    plt.legend()
    plt.grid(True)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    file_to_save = f'{filename}.png'
    plt.savefig(file_to_save, dpi=400, bbox_inches='tight')
    plt.close()

    return bulk_modulus, bulk_modulus_derivative




def calculate_bulk_modulus2(volumes, pressures, filename, dim, start_percent=5, end_percent=95, increment_percent=5):
    """
    Calculates the bulk modulus from pressure-volume data.

    Args:
    volumes (array-like): Array of volume values.
    pressures (array-like): Array of pressure values.
    filename (str): Filename for saving the plot.
    start_percent (int): Starting percentage for segment analysis.
    end_percent (int): Ending percentage for segment analysis.
    increment_percent (int): Increment percentage for segment analysis.

    Returns:
    float: Bulk modulus (K)
    float: Bulk modulus pressure derivative (K')
    """

    def birchmurnaghan_pressure(V, B0, BP, V0):
        eta = (V0 / V)**(1 / 3)
        P = 3 * B0 * (eta**7 - eta**5) * (1 + 0.75 * (BP - 4) * (eta**2 - 1))
        return P


    volumes = np.array(volumes, dtype=float)
    pressures = np.array(pressures, dtype=float)

    if np.all(volumes == volumes[0]):
        print("Error: All volume values are the same. Cannot calculate bulk modulus.")
        return None, None

    def find_best_linear_fit(volumes, pressures, segment_percentages):
        best_mse = float('inf')
        best_segment = 0
        for percent in segment_percentages:
            segment = int(len(volumes) * percent / 100)
            if segment > 5:  # Ensuring enough points for regression
                segment_volumes = volumes[:segment]
                segment_pressures = pressures[:segment]
                slope, intercept, r_value, p_value, std_err = linregress(segment_volumes, segment_pressures)
                mse = np.mean((segment_pressures - (slope * segment_volumes + intercept))**2)
                if mse < best_mse:
                    best_mse = mse
                    best_segment = segment
        return best_segment

    try:
        segment_percentages = range(start_percent, end_percent, increment_percent)
        optimal_segment = find_best_linear_fit(volumes, pressures, segment_percentages)

        linear_volumes = volumes[:optimal_segment] if optimal_segment > 0 else volumes
        linear_pressures = pressures[:optimal_segment] if optimal_segment > 0 else pressures

        if len(linear_volumes) == 0 or len(linear_pressures) == 0:
            raise ValueError("Linear volumes or pressures array is empty.")

        slope, intercept, _, _, _ = linregress(linear_volumes, linear_pressures)
        average_volume = np.mean(linear_volumes)
        bulk_modulus = average_volume * slope

        p0 = [bulk_modulus, 4, average_volume]
        popt, pcov = curve_fit(birchmurnaghan_pressure, volumes, pressures, p0, maxfev=10000)
        B0, BP, V0 = popt
        bulk_modulus_derivative = BP

    except ValueError as e:
        print(f"An error occurred: {e}")
        bulk_modulus = None
        bulk_modulus_derivative = None

    plt.figure(figsize=(10, 6))
    if dim=="2D":
        plt.plot(volumes, pressures, 'o-', label='Pressure-Area Curve')
    else:
        plt.plot(volumes, pressures, 'o-', label='Pressure-Volume Curve')
    if optimal_segment > 0:
        plt.plot(linear_volumes, intercept + slope * linear_volumes, 'r--', label='Linear Fit')

    if dim == "2D":
        plt.xlabel('Area (Å²)', fontsize=20)
        plt.ylabel('Pressure (N/m)', fontsize=20)
        plt.title('Pressure-Area Curve with Linear Fit')
        textstr = f'Stiffness Constant (K): {bulk_modulus:.2f} \nBp (K\'): {bulk_modulus_derivative:.2f}'
    else:
        plt.xlabel('Volume (Å³)', fontsize=20)
        plt.ylabel('Pressure (GPa)', fontsize=20)
        plt.title('Pressure-Volume Curve with Linear Fit')
        textstr = f'Bulk Modulus (K): {bulk_modulus:.2f} GPa\nBp (K\'): {bulk_modulus_derivative:.2f}'
    
    plt.gca().text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=14,
                   verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))

    plt.legend()
    plt.grid(True)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    plt.savefig(f'{filename}.png', dpi=400, bbox_inches='tight')
    plt.close()

    return bulk_modulus, bulk_modulus_derivative


    
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
    elif stress_component_list == "Tensile_biaxial":
        return  (stress_components[0, 0] + stress_components[1, 1]) * conversion_lz/2.
    elif stress_component_list == "xx_yy_xy":
        return  (stress_components[0, 0] + stress_components[1, 1] + stress_components[0, 1])/3. * conversion_lz
    elif stress_component_list == "indent_strength":
        if dim =="2D":
            average_stress = conversion_lz *(stress_components[0, 1] + stress_components[1, 1])/2. 
            #try:
            #    average_stress = conversion_lz * (stress_components[0, 1] + stress_components[1, 1]) / 2.
            #    if average_stress < 0:  # Check if the calculated average stress is negative
            #        raise ValueError("Calculated average stress is negative.")
            #except (IndexError, ValueError):
            #    average_stress = conversion_lz * stress_components[0, 1]
            #average_stress = -conversion_lz*np.sqrt(stress_components[0, 0]**2 - 
            #                          stress_components[0, 0]*stress_components[1, 1] + 
            #                          stress_components[1, 1]**2 + 
            #                          3*stress_components[0, 1]**2)                         
            sigma = np.array([[stress_components[0, 0], stress_components[0, 1]],  # Use eigenvalue method
                              [stress_components[0, 1], stress_components[1, 1]]])

            #average_stress = -conversion_lz*np.linalg.eigvalsh(sigma) 
            #print("DDDDDDDDDDDDDDD", average_stress,aa)
            #average_stress = np.sqrt(average_stress[0]**2 + average_stress[1]**2)
        elif dim =="1D":
            #average_stress = stress_components[0, 1]
            print("Indentation not implemented/debugged for 1D nanotube")
            exit(0)
        elif dim =="3D":
            average_stress =  conversion_lz*(stress_components[0, 2] + stress_components[2, 2])/2. 
           # sigma = np.array([[stress_components[0, 0],stress_components[0, 1], stress_components[2, 0]],  #Use eignenvalue method
           #                   [stress_components[0, 1],stress_components[1, 1], stress_components[1, 2]],
           #                   [stress_components[2, 0],stress_components[1, 2], stress_components[2, 2]]
           #                 ]  )
            
            #average_stress = -conversion_lz*np.linalg.eigvalsh(sigma) 
            #average_stress = np.sqrt(average_stress[0]**2 + average_stress[1]**2 + average_stress[2]**2  )
        return average_stress             
    elif stress_component_list == "ideal_strength":
        if dim == "3D":
            average_stress = stress_components[0, 2] 
        elif dim == "2D":
            average_stress = (stress_components[0, 1] )
        elif dim == "1D":
            average_stress = stress_components[0, 2]

        return average_stress * conversion_lz
    else:
        raise ValueError(f"Unknown stress component: {stress_component_list}")
