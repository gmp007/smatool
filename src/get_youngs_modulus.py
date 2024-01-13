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

def get_youngs_modulusold(stresses, strains):
    num_points = min(5, len(stresses), len(strains))
    
    if num_points == 0:
        raise ValueError("Insufficient data points for computing Young's modulus")

    # Ensure strains and stresses are numpy arrays
    strains_array = np.array(strains[:num_points])  
    stresses_array = np.array(stresses[:num_points]).astype(float)

    mask = (~np.isnan(stresses_array) & ~np.isinf(stresses_array) & 
            ~np.isnan(strains_array) & ~np.isinf(strains_array))

    stresses_array = stresses_array[mask]
    strains_array = strains_array[mask]

    if len(stresses_array) < 2 or len(strains_array) < 2:
        print("Insufficient valid data points after filtering!")
        return None

    A = strains_array[:, np.newaxis]
                    
    #slope, _, _, _ = lstsq(A, stresses_array)
    slope, _, _, _ = np.linalg.lstsq(A, stresses_array, rcond=None)
    elastic_modulus = slope[0]*-0.10 #Conversion factor

    return elastic_modulus
    


def get_youngs_modulus(stresses, strains, window_size=5):
    num_points = min(10, len(stresses), len(strains))
    strains_array = np.array(strains[:num_points])  
    stresses_array = np.array(stresses[:num_points]).astype(float)

    mask = (~np.isnan(stresses_array) & ~np.isinf(stresses_array) & 
            ~np.isnan(strains_array) & ~np.isinf(strains_array))

    stresses = stresses_array[mask]
    strains = strains_array[mask]
    
    
    if len(stresses) != len(strains):
        raise ValueError("Stresses and strains arrays must be of the same length.")

    if len(stresses) < window_size:
        raise ValueError("Insufficient data points for the specified window size.")

    best_r_squared = 0
    best_slope = None

    # Sliding window to find the most linear part of the curve
    for i in range(len(stresses) - window_size + 1):
        window_stresses = stresses[i:i + window_size]
        window_strains = strains[i:i + window_size]

        # Linear regression within the window
        A = np.vstack([window_strains, np.ones(len(window_strains))]).T
        m, c = np.linalg.lstsq(A, window_stresses, rcond=None)[0]

        # Calculate R² (coefficient of determination)
        ss_res = np.sum((window_stresses - (m * window_strains + c)) ** 2)
        ss_tot = np.sum((window_stresses - np.mean(window_stresses)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)

        if r_squared > best_r_squared:
            best_r_squared = r_squared
            best_slope = m*-0.1 # Conversion factor

    if best_slope is not None:
        return best_slope
    else:
        raise ValueError("Could not determine a linear regime in the data.")



