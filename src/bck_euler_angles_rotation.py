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
def slip_systems(hkl,dim,atoms):


    try:
        lattice = atoms.get_cell().array  # Ensure lattice is a numpy array
        positions = atoms.get_scaled_positions()  # Scaled positions
        numbers = atoms.get_atomic_numbers()  # Atomic numbers
        cell = (lattice, positions, numbers)
        spg = spglib.get_spacegroup(cell, symprec=0.1)
    except TypeError:
        lattice = atoms.get_cell().array
        positions = atoms.get_scaled_positions()
        numbers = atoms.get_atomic_numbers()
        cell = (lattice, positions, numbers)
        spg = spglib.get_spacegroup(cell, symprec=0.1)

    # Extract the space group number
    spg_number = int(spg.split()[1].strip('()'))
    
    phi = phi1 = phi2 = None 

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
        phi1, phi, phi2 = 0.,0.,0.
        if len(hkl) > 2:
            print("Warning: Only 'h' and 'k' values are used for 2D crystal systems. l will be ignored.")
        #h, k = hkl
        h, k = hkl[:2]
    
        if crystal_system == 'hexagonal':

            norm = np.sqrt(h**2 + k**2)
            phi1 = np.arctan2(np.sqrt(3) * k / 2, h) * 180 / np.pi

        if crystal_system == 'tetragonal':

            a_comp = np.sqrt(h**2 + k**2)
            phi1 = np.arctan2(k, h) * 180 / np.pi

            #phi2 = phi1 + 90 * np.round((phi1 - 45) / 90)
        if crystal_system == 'trigonal':

            a_comp = (2 * k - h) / np.sqrt(3)
            #c_comp = l
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

    lattice_vectors = np.array([list(map(float, line.split())) for line in lines[2:5]])

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
    strain_alpha, strain_beta, strain_gamma = slip_systems(strain_direction, dim, atoms)
    print(f"Euler Angles for the Strain Direction: alpha:{strain_alpha:.3f}, beta:{strain_beta:.3f}, gamma:{strain_gamma:.3f}" )
    if slipon:
        slip_alpha, slip_beta, slip_gamma = slip_systems(slip_direction, dim, atoms)
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
    alpha = float(alpha) * pi / 180  

    # Rotate around z-axis in 2D
    for i in range(len(matrix)):
        x, y = matrix[i][0], matrix[i][1]
        matrix[i][0] = x * math.cos(alpha) - y * math.sin(alpha)
        matrix[i][1] = x * math.sin(alpha) + y * math.cos(alpha)

    return matrix
    
def rotate_crystal_structure_2D(poscar_file, strain_direction, slip_direction, dim, atoms,slipon):
    strain_phi,_,_ = slip_systems(strain_direction, dim, atoms)
    print(f"Euler Angle for Strain Direction. For 2D only z-rotation needed: gamma:{strain_phi:.3f}" )
    if slipon:
        slip_phi,_,_ =   slip_systems(slip_direction, dim, atoms)
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
    

def generate_slip_systems_2d_hexagonalold():
    """
    Generate slip systems for a 2D hexagonal crystal system (like graphene).
    Returns a list of slip directions.
    """
    # Hexagonal lattice vectors
    a1 = np.array([1, 0])  # along x-axis
    a2 = np.array([np.cos(np.pi/3), np.sin(np.pi/3)])  # 60 degrees from x-axis

    # Define slip directions for a 2D hexagonal system
    # Three primary slip directions based on hexagonal symmetry
    directions = []
    directions = [
        a1,  # along the first lattice vector
        -0.5 * a1 + np.sqrt(3)/2 * a2,  # 120 degrees from a1
        -0.5 * a1 - np.sqrt(3)/2 * a2   # -120 degrees from a1
    ]

    return directions

def generate_slip_systems_2d_hexagonal():
    """
    We assume hexagonal for 2D materials
    Generate slip systems for a 2D hexagonal crystal system.
    :return: List of tuples containing normal vectors and slip directions in 2D
    """
    # Define slip directions for 2D hexagonal system
    # Assuming two primary slip directions for simplicity
    directions = [(1, 0), (0, 1), (np.cos(np.pi/3), np.sin(np.pi/3))]  # Example directions

    # In 2D, the normal to the slip plane is not explicitly defined as in 3D
    # For simplicity, we can assume the plane normal is perpendicular to the slip direction
    slip_systems = []
    for direction in directions:
        normal_vector = (-direction[1], direction[0])  # Rotate by 90 degrees to get normal
        slip_systems.append((normal_vector, direction))

    return slip_systems



def get_cubic_crystal_type(spg_number):
    """Determine the crystal type for cubic space groups."""
    if 195 <= spg_number <= 199:
        return "sc"  # Simple Cubic
    elif 200 <= spg_number <= 206:
        return "bcc"  # Body-Centered Cubic
    elif 207 <= spg_number <= 230:
        return "fcc"  # Face-Centered Cubic
    else:
        raise ValueError("Unsupported cubic space group: " + str(spg_number))

        
        
def get_cubic_crystal_type(spg_number):
    """
    Determine the cubic crystal type (FCC, BCC, SC) based on the space group number.

    Args:
        spg_number: Integer representing the space group number.

    Returns:
        String indicating the cubic crystal type (FCC, BCC, SC).
    """

    cubic_types = {
        195: "sc",   # P23
        196: "sc",   # F23
        197: "sc",   # I23
        198: "sc",   # P213
        199: "sc",   # I213
        200: "bcc",  # Pm-3
        201: "bcc",  # Pn-3
        202: "bcc",  # Fm-3
        203: "bcc",  # Fd-3
        204: "bcc",  # Im-3
        205: "bcc",  # Pa-3
        206: "bcc",  # Ia-3
        207: "fcc",  # Pm-3m
        208: "fcc",  # Pn-3n
        209: "fcc",  # Fm-3m
        210: "fcc",  # Fd-3m
        211: "fcc",  # Im-3m
        212: "fcc",  # Pa-3
        213: "fcc",  # Ia-3
        214: "fcc",  # P-43m
        215: "fcc",  # I-43m
        216: "fcc",  # P-43n
        217: "fcc",  # F-43m
        218: "fcc",  # Fd-3m
        219: "fcc",  # Im-3m
        220: "fcc",  # Pa-3
        221: "fcc",  # Ia-3
        222: "fcc",  # Pm-3n
        223: "fcc",  # Pn-3n
        224: "fcc",  # Pm-3m
        225: "fcc",  # Pn-3n
        226: "fcc",  # Pm-3m
        227: "fcc",  # Pn-3n
        228: "fcc",  # Pm-3m
        229: "fcc",  # Pn-3n
        230: "fcc"   # Pm-3m
    }

    return cubic_types.get(spg_number, "unknown")



def generate_slip_systems(atoms):
    """
    Generate slip systems for specified crystal systems.
    :param atoms: ASE Atoms object
    :return: List of tuples containing normal vectors and slip directions
    """
    # Determine the space group
    try:
        lattice = atoms.get_cell().array  # Ensure lattice is a numpy array
        positions = atoms.get_scaled_positions()  # Scaled positions
        numbers = atoms.get_atomic_numbers()  # Atomic numbers
        cell = (lattice, positions, numbers)
        spg = spglib.get_spacegroup(cell, symprec=0.1)
    except TypeError:
        lattice = atoms.get_cell().array
        positions = atoms.get_scaled_positions()
        numbers = atoms.get_atomic_numbers()
        cell = (lattice, positions, numbers)
        spg = spglib.get_spacegroup(cell, symprec=0.1)

    # Extract the space group number
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
        crystal_type = get_cubic_crystal_type(spg_number)
        if crystal_type == "fcc":
          # FCC: {111} planes with  directions
          planes = [(1, 1, 1), (-1, 1, 1), (1, -1, 1), (1, 1, -1)]
          directions = [(1, 1, 0), (-1, 1, 0), (1, -1, 0), (1, 0, 1), (-1, 0, 1), (0, 1, 1), (0, -1, 1)]
        elif crystal_type == "bcc":
          # BCC: {110}, {123}, and {112} planes with  directions
          planes = [(1, 1, 0), (-1, 1, 0), (1, -1, 0), (0, 1, 1), (0, -1, 1), (1, 0, 1), (-1, 0, 1),
                    (1, 2, 3), (-1, 2, 3), (1, -2, 3), (1, 2, -3), (-1, -2, 3), (1, -2, -3), (-1, 2, -3),
                    (1, 1, 2), (-1, 1, 2), (1, -1, 2), (1, 1, -2), (-1, -1, 2), (1, -1, -2), (-1, 1, -2)]
          directions = [(1, 1, 1), (-1, 1, 1), (1, -1, 1), (-1, -1, 1)]
        elif crystal_type == "sc":
          # SC: {100}, {010}, and {001} planes with  directions
          planes = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
          directions = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        else:
          raise ValueError("Unsupported cubic crystal type: " + crystal_type)


        
    elif crystal_system == 'hexagonal':
        # Define basal, prismatic, and pyramidal planes and <11-20> directions for hexagonal crystals
        basal_planes = [(0, 0, 1), (0, 0, -1)]
        prismatic_planes = [(1, 0, -1), (-1, 0, 1), (0, 1, -1), (0, -1, 1)]
        pyramidal_planes = [(1, 0, -1), (-1, 0, 1), (0, 1, -1), (0, -1, 1), (1, -1, 0), (-1, 1, 0)]
        directions = [(1, 1, -2), (-1, -1, 2)]
        planes = basal_planes + prismatic_planes + pyramidal_planes
        #slip_systems = [(plane, direction) for plane in planes for direction in directions]
        
    elif crystal_system == 'triclinic':
        # Triclinic crystals have no symmetry constraints, so defining slip systems is more complex and less standardized
        index_range = [-1, 0, 1]
        # Generate all possible combinations of indices for planes and directions
        planes = list(itertools.product(index_range, repeat=3))
        directions = list(itertools.product(index_range, repeat=3))

        # Remove the (0, 0, 0) vector as it's not valid for planes or directions
        planes = [p for p in planes if p != (0, 0, 0)]
        directions = [d for d in directions if d != (0, 0, 0)]
        #slip_systems = [(plane, direction) for plane in planes for direction in directions]
        
    elif crystal_system == 'tetragonal':
        planes = [(1, 1, 0), (1, -1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 2), (1, -1, 2), (0, 0, 1)]
        directions = [(1, 1, 1), (-1, 1, 1), (0, 0, 1)]
        #slip_systems = [(plane, direction) for plane in planes for direction in directions]
        
    elif crystal_system == 'trigonal':
        basal_planes = [(0, 0, 1), (0, 0, -1)]
        prismatic_planes = [(1, 0, -1), (-1, 0, 1), (0, 1, -1), (0, -1, 1)]
        pyramidal_planes = [(1, 0, -1), (-1, 0, 1), (0, 1, -1), (0, -1, 1), (1, -1, 0), (-1, 1, 0)]
        directions = [(1, 1, -2), (-1, -1, 2)]
        planes = basal_planes + prismatic_planes + pyramidal_planes
        #slip_systems = [(plane, direction) for plane in planes for direction in directions]
            
    elif crystal_system == 'orthorhombic':
        planes = [(1, 1, 0), (0, 1, 1), (1, 0, 1), (0, 0, 1)]
        directions = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        #slip_systems = [(plane, direction) for plane in planes for direction in directions]
        
    elif crystal_system == 'monoclinic':

        planes = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1)]
        directions = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        #slip_systems = [(plane, direction) for plane in planes for direction in directions]
    else:
        raise ValueError("Unsupported crystal system.")

    # Generate slip systems
    slip_systems = []
    for plane in planes:
        for direction in directions:
            if np.dot(plane, direction) == 0:  # Ensure the direction is within the plane
                slip_systems.append((plane, direction))
    # Remove duplicates
    slip_systems = list(set(slip_systems))
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
