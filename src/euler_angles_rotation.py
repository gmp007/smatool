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
import spglib
import math
import itertools
from read_write import read_options_from_input

options = read_options_from_input()
slipon = options.get("slipon",False)
def calculate_euler_angles(hkl,dim,atoms):

    # Determine the space group
    spg = spglib.get_spacegroup(atoms,symprec=0.1)
    spg_number = int(spg.split()[1].strip('()'))
    phi = phi2 = None 

    if 1 <= spg_number <= 2:
        crystal_system = "triclinic"
    elif 3 <= spg_number <= 15:
        crystal_system = "monoclinic"
    elif 16 <= spg_number <= 74:
        crystal_system = "orthorhombic"
    elif 75 <= spg_number <= 142:
        crystal_system = "tetragonal"
    elif 143 <= spg_number <= 167:
        crystal_system = "trigonal"
    elif 168 <= spg_number <= 194:
        crystal_system = "hexagonal"
    elif 195 <= spg_number <= 230:
        crystal_system = "cubic"
    else:
        crystal_system = "unknown"  

    if dim =="3D" or dim == "1D":

        h, k, l = hkl

        if crystal_system == 'cubic':

            norm = np.sqrt(h**2 + k**2 + l**2)
            a_cos = h / norm
            b_cos = k / norm
            c_cos = l / norm


            x = a_cos * b_cos + c_cos * (-np.sqrt(1 - a_cos**2))
            y = b_cos * c_cos + a_cos * (-np.sqrt(1 - b_cos**2))

            phi1 = np.arctan2(y, x) * 180 / np.pi
            phi = np.arccos(c_cos) * 180 / np.pi
            phi2 = np.arctan2(a_cos * (-np.sqrt(1 - c_cos**2)) - b_cos * c_cos,
                              np.sqrt(1 - a_cos**2 - b_cos**2)) * 180 / np.pi

        elif crystal_system == 'hexagonal':

            a_comp = np.sqrt(3) * k / 2
            c_comp = l

            phi1 = np.arctan2(a_comp, h) * 180 / np.pi

            phi = np.arccos(c_comp / np.sqrt(h**2 + k**2 + l**2)) * 180 / np.pi

            phi2 = phi1 + 120 * np.round((phi1 - 30) / 120)

        elif crystal_system == 'trigonal':

            a_comp = (2 * k - h) / np.sqrt(3)
            c_comp = l


            phi1 = np.arctan2(a_comp, np.sqrt(3) * h) * 180 / np.pi
            phi = np.arccos(c_comp / np.sqrt(h**2 + k**2 + l**2)) * 180 / np.pi
            phi2 = phi1 + 120 * np.round((phi1 - 30) / 120)

        elif crystal_system == 'tetragonal':

            a_comp = np.sqrt(h**2 + k**2)
            c_comp = l

            phi1 = np.arctan2(k, h) * 180 / np.pi
            phi = np.arccos(c_comp / np.sqrt(a_comp**2 + c_comp**2)) * 180 / np.pi
            phi2 = phi1 + 180 * np.round((phi1 - 90) / 180)

        elif crystal_system == 'orthorhombic':

            norm = np.sqrt(h**2 + k**2 + l**2)

            a_cos = h / norm
            b_cos = k / norm
            c_cos = l / norm

            x = a_cos * np.sqrt(1 - c_cos**2)
            y = b_cos * np.sqrt(1 - c_cos**2)

            phi1 = np.arctan2(y, x) * 180 / np.pi
            phi = np.arccos(c_cos) * 180 / np.pi
            phi2 = np.arctan2(a_cos * c_cos, b_cos * c_cos) * 180 / np.pi

        elif crystal_system == 'monoclinic':
            norm = np.sqrt(h**2 + k**2 + l**2)
            a_cos = h / norm
            b_cos = k / norm
            c_cos = l / norm

            beta = np.arccos(b_cos) * 180 / np.pi

            x_bar = a_cos * np.cos(beta) + c_cos * np.sin(beta)
            y_bar = k / norm
            z_bar = -a_cos * np.sin(beta) + c_cos * np.cos(beta)

            phi1 = np.arctan2(y_bar, x_bar) * 180 / np.pi
            phi2 = np.arctan2(np.sqrt(1 - y_bar**2), y_bar) * 180 / np.pi
            phi = np.arccos(z_bar) * 180 / np.pi

        elif crystal_system == 'triclinic':

            norm = np.sqrt(h**2 + k**2 + l**2)
            a_cos = h / norm
            b_cos = k / norm
            c_cos = l / norm

            x = a_cos * b_cos * c_cos - c_cos * np.sqrt(1 - a_cos**2) * np.sqrt(1 - b_cos**2)
            y = a_cos * np.sqrt(1 - c_cos**2) + b_cos * c_cos * np.sqrt(1 - a_cos**2)
            z = -b_cos * np.sqrt(1 - c_cos**2) + a_cos * c_cos * np.sqrt(1 - b_cos**2)

            phi1 = np.arctan2(y, x) * 180 / np.pi
            phi = np.arccos(z) * 180 / np.pi
            phi2 = np.arctan2(a_cos * c_cos * np.sqrt(1 - b_cos**2) - b_cos * np.sqrt(1 - a_cos**2),
                                b_cos * c_cos * np.sqrt(1 - a_cos**2) + a_cos * np.sqrt(1 - b_cos**2)) * 180 / np.pi

        else:
            raise ValueError("Unsupported crystal system.")

    elif dim =="2D":
        
        if len(hkl) > 2:
            print("Warning: Only 'h' and 'k' values are used for 2D crystal systems. Please ensure this is correct.")
        h, k = hkl
    
        if crystal_system == 'hexagonal':

            norm = np.sqrt(h**2 + k**2)
            a_cos = h / norm
            phi1 = np.arctan2(np.sqrt(3) * k / 2, h) * 180 / np.pi

        if crystal_system == 'tetragonal':

            a_comp = np.sqrt(h**2 + k**2)

            phi1 = np.arctan2(k, h) * 180 / np.pi

            #phi2 = phi1 + 90 * np.round((phi1 - 45) / 90)
        if crystal_system == 'trigonal':

            a_comp = (2 * k - h) / np.sqrt(3)
            c_comp = l

            # Calculate in-plane angle phi1
            phi1 = np.arctan2(a_comp, np.sqrt(3) * h) * 180 / np.pi

        else:

            a_comp = np.sqrt(h**2 + k**2)
            phi1 = np.arctan2(k, h) * 180 / np.pi

    else:
        raise ValueError("Invalid dimensionality. Please specify '1D', '2D', or '3D'.")


    return phi1, phi, phi2
    



def apply_rotation(matrix, alpha, beta, gamma):
    """
    Applies a rotation to a matrix.

    Args:
      matrix: A 2D NumPy array representing the points to rotate.
      alpha: Rotation angle around the Z-axis in degrees.
      beta: Rotation angle around the Y-axis in degrees.
      gamma: Rotation angle around the X-axis in degrees.

    Returns:
      The rotated matrix.
    """
    # Convert angles from degrees to radians
    alpha, beta, gamma = np.radians([alpha, beta, gamma])

    # Rotation matrices
    Rz = np.array([[math.cos(alpha), -math.sin(alpha), 0],
                   [math.sin(alpha), math.cos(alpha), 0],
                   [0, 0, 1]])
    Ry = np.array([[math.cos(beta), 0, math.sin(beta)],
                   [0, 1, 0],
                   [-math.sin(beta), 0, math.cos(beta)]])
    Rx = np.array([[1, 0, 0],
                   [0, math.cos(gamma), -math.sin(gamma)],
                   [0, math.sin(gamma), math.cos(gamma)]])

    # Apply rotations
    rotated_matrix = np.dot(matrix, Rz)
    rotated_matrix = np.dot(rotated_matrix, Ry)
    rotated_matrix = np.dot(rotated_matrix, Rx)

    return rotated_matrix



def rotate(pos, alpha, beta, gamma):
    with open(pos) as f:
        lines = f.readlines()

    # Extract lattice vectors
    lattice_vectors = np.array([list(map(float, line.split())) for line in lines[2:5]])

    # Convert angles from degrees to radians
    alpha, beta, gamma = np.radians([alpha, beta, gamma])

    # Rotation matrices
    Rz = np.array([[math.cos(alpha), -math.sin(alpha), 0],
                   [math.sin(alpha), math.cos(alpha), 0],
                   [0, 0, 1]])
    Ry = np.array([[math.cos(beta), 0, math.sin(beta)],
                   [0, 1, 0],
                   [-math.sin(beta), 0, math.cos(beta)]])
    Rx = np.array([[1, 0, 0],
                   [0, math.cos(gamma), -math.sin(gamma)],
                   [0, math.sin(gamma), math.cos(gamma)]])

    # Apply rotations
    rotated_vectors = np.dot(lattice_vectors, Rz)
    rotated_vectors = np.dot(rotated_vectors, Ry)
    rotated_vectors = np.dot(rotated_vectors, Rx)
    print(rotated_vectors)
    # Write to file
    with open("POSCAR_rotated", 'w') as f:
        f.writelines(lines[:2])
        for vec in rotated_vectors:
            f.write(' '.join(map(str, vec)) + '\n')
        f.writelines(lines[5:])
        
        
    
def rotate_crystal_structure(poscar_file, strain_direction, slip_direction, dim, atoms,slipon):
    strain_alpha, strain_beta, strain_gamma = calculate_euler_angles(strain_direction, dim, atoms)
    print(f"Euler Angles for the Strain Direction: alpha:{strain_alpha:.3f}, beta:{strain_beta:.3f}, gamma:{strain_gamma:.3f}" )
    if slipon:
        slip_alpha, slip_beta, slip_gamma = calculate_euler_angles(slip_direction, dim, atoms)
        print(f"Euler Angles for the Slip Direction: alpha:{slip_alpha:.3f}, beta:{slip_beta:.3f}, gamma:{slip_gamma:.3f}" )

    with open(poscar_file, 'r') as file:
        lines = file.readlines()

    lattice_vectors = np.array([list(map(float, line.split())) for line in lines[2:5]])

    lattice_vectors = apply_rotation(lattice_vectors, strain_alpha, strain_beta, strain_gamma)
    if slipon:
        lattice_vectors = apply_rotation(lattice_vectors, slip_alpha, slip_beta, slip_gamma)

    
    with open("POSCAR_rotated", 'w') as file:
        file.writelines(lines[:2])
        for vec in lattice_vectors:
            file.write(' '.join(map(str, vec)) + '\n')
        file.writelines(lines[5:])

    print("Structure rotated and saved in VASP format as POSCAR_rotated.")



def apply_rotation_2D(matrix, alpha):
    pi = np.pi
    alpha = float(alpha) * pi / 180  # Convert to radians

    # Rotate around z-axis in 2D
    for i in range(len(matrix)):
        x, y = matrix[i][0], matrix[i][1]
        matrix[i][0] = x * math.cos(alpha) - y * math.sin(alpha)
        matrix[i][1] = x * math.sin(alpha) + y * math.cos(alpha)

    return matrix
    
def rotate_crystal_structure_2D(poscar_file, strain_direction, slip_direction, dim, atoms,slipon):
    strain_phi,_,_ = calculate_euler_angles(strain_direction, dim, atoms)
    print(f"Euler Angle for Strain Direction. For 2D only z-rotation needed: gamma:{strain_phi:.3f}" )
    if slipon:
        slip_phi,_,_ =   calculate_euler_angles(slip_direction, dim, atoms)
        print(f"Euler Angle for Slip Direction For 2D only z-rotation needed: gamma:{slip_phi:.3f}" )

    
    with open(poscar_file, 'r') as file:
        lines = file.readlines()

    lattice_vectors = np.array([list(map(float, line.split()))[:2] for line in lines[2:4]])
    lattice_vectors = apply_rotation_2D(lattice_vectors, strain_phi)
    
    if slipon:
        lattice_vectors = apply_rotation_2D(lattice_vectors, slip_phi)

    with open("POSCAR_rotated", 'w') as f:
        for line in lines[:2]:
            f.write(line)
        for vec in lattice_vectors:
            f.write(f"{vec[0]:.6f} {vec[1]:.6f} 0.0\n")
        for line in lines[4:]:
            f.write(line)
            
    print("Structure rotated and saved in VASP format as POSCAR_rotated.")
    

def generate_slip_systems_2d_hexagonal():
    """
    We assume hexagonal for 2D materials
    Generate slip systems for a 2D hexagonal crystal system.
    :return: List of tuples containing normal vectors and slip directions in 2D
    """
    # Define slip directions for 2D hexagonal system
    # Assuming two primary slip directions for simplicity
    directions = [(1, 0), (0, 1)]  # Example directions

    # In 2D, the normal to the slip plane is not explicitly defined as in 3D
    # For simplicity, we can assume the plane normal is perpendicular to the slip direction
    slip_systems = []
    for direction in directions:
        normal_vector = (-direction[1], direction[0])  # Rotate by 90 degrees to get normal
        slip_systems.append((normal_vector, direction))

    return slip_systems



def generate_slip_systems(atoms):
    """
    Generate slip systems for specified crystal systems.
    :param crystal_system: Type of crystal system ('cubic' or 'hexagonal')
    :return: List of tuples containing normal vectors and slip directions
    """
    # Determine the space group
    spg = spglib.get_spacegroup(atoms,symprec=0.1)
    spg_number = int(spg.split()[1].strip('()'))
    if 1 <= spg_number <= 2:
        crystal_system = "triclinic"
    elif 3 <= spg_number <= 15:
        crystal_system = "monoclinic"
    elif 16 <= spg_number <= 74:
        crystal_system = "orthorhombic"
    elif 75 <= spg_number <= 142:
        crystal_system = "tetragonal"
    elif 143 <= spg_number <= 167:
        crystal_system = "trigonal"
    elif 168 <= spg_number <= 194:
        crystal_system = "hexagonal"
    elif 195 <= spg_number <= 230:
        crystal_system = "cubic"
    else:
        crystal_system = "unknown"   
    if crystal_system == 'cubic':
        # Define {111} planes and <110> directions for cubic crystals
        planes = [(1, 1, 1), (-1, 1, 1), (1, -1, 1), (1, 1, -1)]
        directions = [(1, 1, 0), (-1, 1, 0), (1, -1, 0), (1, 0, 1), (-1, 0, 1), (0, 1, 1), (0, -1, 1)]

    elif crystal_system == 'hexagonal':
        # Define basal, prismatic, and pyramidal planes and <11-20> directions for hexagonal crystals
        basal_planes = [(0, 0, 1), (0, 0, -1)]
        prismatic_planes = [(1, 0, -1), (-1, 0, 1), (0, 1, -1), (0, -1, 1)]
        pyramidal_planes = [(1, 0, -1), (-1, 0, 1), (0, 1, -1), (0, -1, 1), (1, -1, 0), (-1, 1, 0)]
        directions = [(1, 1, -2), (-1, -1, 2)]
        planes = basal_planes + prismatic_planes + pyramidal_planes

    elif crystal_system == 'triclinic':
        # Triclinic crystals have no symmetry constraints, so defining slip systems is more complex and less standardized
        # Define a range for indices to generate planes and directions
        index_range = [-1, 0, 1]
        # Generate all possible combinations of indices for planes and directions
        planes = list(itertools.product(index_range, repeat=3))
        directions = list(itertools.product(index_range, repeat=3))

        # Remove the (0, 0, 0) vector as it's not valid for planes or directions
        planes = [p for p in planes if p != (0, 0, 0)]
        directions = [d for d in directions if d != (0, 0, 0)]

    elif crystal_system == 'tetragonal':
        planes = [(1, 1, 0), (1, -1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 2), (1, -1, 2), (0, 0, 1)]
        directions = [(1, 1, 1), (-1, 1, 1), (0, 0, 1)]

    elif crystal_system == 'trigonal':
        basal_planes = [(0, 0, 1), (0, 0, -1)]
        prismatic_planes = [(1, 0, -1), (-1, 0, 1), (0, 1, -1), (0, -1, 1)]
        pyramidal_planes = [(1, 0, -1), (-1, 0, 1), (0, 1, -1), (0, -1, 1), (1, -1, 0), (-1, 1, 0)]
        directions = [(1, 1, -2), (-1, -1, 2)]
        planes = basal_planes + prismatic_planes + pyramidal_planes
    
    elif crystal_system == 'orthorhombic':
        planes = [(1, 1, 0), (0, 1, 1), (1, 0, 1), (0, 0, 1)]
        directions = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

    elif crystal_system == 'monoclinic':

        planes = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1)]
        directions = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

    else:
        raise ValueError("Unsupported crystal system.")

    # Generate slip systems
    slip_systems = []
    for plane in planes:
        for direction in directions:
            if np.dot(plane, direction) == 0:  # Ensure the direction is within the plane
                slip_systems.append((plane, direction))

    return slip_systems,crystal_system



def calculate_schmid_factor(normal_vector, slip_direction, stress_direction):
    """
    Calculate the Schmid factor for a given slip system.

    :param normal_vector: Normal vector to the slip plane (as a list or tuple of 3 elements)
    :param slip_direction: Slip direction (as a list or tuple of 3 elements)
    :param stress_direction: Direction of the applied stress (as a list or tuple of 3 elements)
    :return: Schmid factor for the given slip system
    """
 #   if len(stress_direction) != 3:
 #       print("Error: Strain direction must have 3 elements.")
 #       return None
    # Normalize the vectors
    normal_vector = np.array(normal_vector) / np.linalg.norm(normal_vector)
    slip_direction = np.array(slip_direction) / np.linalg.norm(slip_direction)
    stress_direction = np.array(stress_direction) / np.linalg.norm(stress_direction)

    # Calculate the Schmid factor
    return np.dot(normal_vector, stress_direction) * np.dot(slip_direction, stress_direction)
