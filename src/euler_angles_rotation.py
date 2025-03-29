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
slipon = options.get("slipon", False)

def slip_systems(hkl, dim, atoms):
    try:
        lattice = atoms.get_cell()  # Lattice vectors
        positions = atoms.get_scaled_positions()  # Scaled positions
        numbers = atoms.get_atomic_numbers()  # Atomic numbers
        cell = (lattice, positions, numbers)
        spg = spglib.get_spacegroup(cell, symprec=0.1)
    except Exception as e:
        raise RuntimeError(f"Error determining space group: {e}")

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

    if dim == "3D" or dim == "1D":
        h, k, l = hkl

        if crystal_system == 'cubic':
            norm = np.sqrt(h**2 + k**2 + l**2)
            a_cos = h / norm
            b_cos = k / norm
            c_cos = l / norm

            # Euler angles for cubic crystals
            phi1 = np.arctan2(b_cos, a_cos) * 180 / np.pi
            phi = np.arccos(c_cos) * 180 / np.pi
            phi2 = 0.0  # For cubic crystals, the third angle can be set to zero

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
            y_bar = b_cos
            z_bar = -a_cos * np.sin(beta) + c_cos * np.cos(beta)

            phi1 = np.arctan2(y_bar, x_bar) * 180 / np.pi
            phi = np.arccos(z_bar) * 180 / np.pi
            phi2 = 0.0  # Simplification

        elif crystal_system == 'triclinic':
            # For triclinic systems, Euler angles can be complex to calculate
            raise NotImplementedError("Euler angle calculation for triclinic systems is not implemented.")
        else:
            raise ValueError("Unsupported crystal system.")

    elif dim == "2D":
        phi1 = phi = phi2 = 0.0
        h, k = hkl[:2]

        if crystal_system == 'hexagonal':
            norm = np.sqrt(h**2 + k**2)
            phi1 = np.arctan2(np.sqrt(3) * k / 2, h) * 180 / np.pi

        elif crystal_system == 'tetragonal' or crystal_system == 'orthorhombic':
            phi1 = np.arctan2(k, h) * 180 / np.pi

        elif crystal_system == 'trigonal':
            a_comp = (2 * k - h) / np.sqrt(3)
            phi1 = np.arctan2(a_comp, np.sqrt(3) * h) * 180 / np.pi
        else:
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
    alpha_rad, beta_rad, gamma_rad = np.radians([alpha, beta, gamma])

    # Rotation matrices
    Rz = np.array([[math.cos(alpha_rad), -math.sin(alpha_rad), 0],
                   [math.sin(alpha_rad), math.cos(alpha_rad), 0],
                   [0, 0, 1]])
    Ry = np.array([[math.cos(beta_rad), 0, math.sin(beta_rad)],
                   [0, 1, 0],
                   [-math.sin(beta_rad), 0, math.cos(beta_rad)]])
    Rx = np.array([[1, 0, 0],
                   [0, math.cos(gamma_rad), -math.sin(gamma_rad)],
                   [0, math.sin(gamma_rad), math.cos(gamma_rad)]])

    # Combined rotation matrix
    R = Rz @ Ry @ Rx

    # Apply rotation
    rotated_matrix = matrix @ R.T  # Transpose R to apply rotation correctly

    return rotated_matrix

def rotate_crystal_structure(poscar_file, strain_direction, slip_direction, dim, atoms, slipon):
    strain_alpha, strain_beta, strain_gamma = slip_systems(strain_direction, dim, atoms)
    print(f"Euler Angles for the Strain Direction: alpha: {strain_alpha:.3f}, beta: {strain_beta:.3f}, gamma: {strain_gamma:.3f}")
    if slipon:
        slip_alpha, slip_beta, slip_gamma = slip_systems(slip_direction, dim, atoms)
        print(f"Euler Angles for the Slip Direction: alpha: {slip_alpha:.3f}, beta: {slip_beta:.3f}, gamma: {slip_gamma:.3f}")

    with open(poscar_file, 'r') as file:
        lines = file.readlines()

    # Extract lattice vectors from POSCAR
    lattice_vectors = np.array([list(map(float, line.strip().split())) for line in lines[2:5]])

    # Ensure lattice vectors are 3D
    if lattice_vectors.shape != (3, 3):
        raise ValueError("Lattice vectors should be a 3x3 matrix.")

    # Apply rotations
    lattice_vectors = apply_rotation(lattice_vectors, strain_alpha, strain_beta, strain_gamma)
    if slipon:
        lattice_vectors = apply_rotation(lattice_vectors, slip_alpha, slip_beta, slip_gamma)

    # Write the rotated structure to a new POSCAR file
    with open("POSCAR_rotated", 'w') as file:
        file.writelines(lines[:2])
        for vec in lattice_vectors:
            file.write('  {:.16f}  {:.16f}  {:.16f}\n'.format(*vec))
        file.writelines(lines[5:])

    print("Structure rotated and saved in VASP format as POSCAR_rotated.")

def apply_rotation_2D(matrix, alpha):
    """
    Applies a rotation around the z-axis in 2D.
    """
    alpha_rad = np.radians(alpha)
    Rz = np.array([[math.cos(alpha_rad), -math.sin(alpha_rad)],
                   [math.sin(alpha_rad), math.cos(alpha_rad)]])

    # Apply rotation
    rotated_matrix = matrix @ Rz.T

    return rotated_matrix

def rotate_crystal_structure_2D(poscar_file, strain_direction, slip_direction, dim, atoms, slipon):
    strain_phi, _, _ = slip_systems(strain_direction, dim, atoms)
    print(f"Euler Angle for Strain Direction (2D rotation): gamma: {strain_phi:.3f}")
    if slipon:
        slip_phi, _, _ = slip_systems(slip_direction, dim, atoms)
        print(f"Euler Angle for Slip Direction (2D rotation): gamma: {slip_phi:.3f}")

    with open(poscar_file, 'r') as file:
        lines = file.readlines()

    # Extract 2D lattice vectors (first two components)
    lattice_vectors = np.array([list(map(float, line.strip().split()))[:2] for line in lines[2:4]])

    # Apply rotations
    lattice_vectors = apply_rotation_2D(lattice_vectors, strain_phi)
    if slipon:
        lattice_vectors = apply_rotation_2D(lattice_vectors, slip_phi)

    # Write the rotated structure to a new POSCAR file
    with open("POSCAR_rotated", 'w') as f:
        f.writelines(lines[:2])
        for vec in lattice_vectors:
            f.write(f"  {vec[0]:.16f}  {vec[1]:.16f}  0.0000000000000000\n")
        f.writelines(lines[4:])

    print("Structure rotated and saved in VASP format as POSCAR_rotated.")

def get_cubic_crystal_type(spg_number):
    """
    Determine the cubic crystal type (FCC, BCC, SC) based on the space group number.
    """
    if 195 <= spg_number <= 199:
        return "sc"  # Simple Cubic
    elif 200 <= spg_number <= 206:
        return "bcc"  # Body-Centered Cubic
    elif 207 <= spg_number <= 230:
        return "fcc"  # Face-Centered Cubic
    else:
        return "unknown"



def generate_slip_systems_2d_hexagonal():
    """
    Generate slip systems for a 2D hexagonal crystal system.
    :return: List of tuples containing normal vectors and slip directions in 2D
    """
    # Define 2D hexagonal lattice vectors
    a1 = np.array([1, 0])  # Along x-axis
    a2 = np.array([0.5, np.sqrt(3)/2])  # 60 degrees from x-axis

    # Slip directions are along the lattice vectors and their negatives
    directions = [
        a1,             # Along a1
        a2,             # Along a2
        -a1,            # Opposite to a1
        -a2,            # Opposite to a2
        a1 - a2,        # Along a3 = a1 - a2
        -(a1 - a2)      # Opposite to a3
    ]

    # For 2D, the normal to the slip plane is perpendicular to the slip direction
    slip_systems = []
    for direction in directions:
        # Normalize the slip direction
        direction_norm = direction / np.linalg.norm(direction)
        # Calculate the normal vector by rotating the slip direction by 90 degrees
        normal_vector = np.array([-direction_norm[1], direction_norm[0]])
        slip_systems.append((normal_vector, direction_norm))

    return slip_systems


def generate_slip_systems(atoms):
    """
    Generate slip systems for specified crystal systems.
    :param atoms: ASE Atoms object
    :return: List of tuples containing normal vectors and slip directions
    """
    # Determine the space group
    try:
        lattice = atoms.get_cell()
        positions = atoms.get_scaled_positions()
        numbers = atoms.get_atomic_numbers()
        cell = (lattice, positions, numbers)
        spg = spglib.get_spacegroup(cell, symprec=0.1)
    except Exception as e:
        raise RuntimeError(f"Error determining space group: {e}")

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

    planes = []
    directions = []

    if crystal_system == 'cubic':
        crystal_type = get_cubic_crystal_type(spg_number)
        if crystal_type == "fcc":
            # FCC: {111} planes with <110> directions
            planes = [(1, 1, 1), (-1, 1, 1), (1, -1, 1), (1, 1, -1)]
            directions = [(1, 1, 0), (-1, 1, 0), (1, -1, 0), (1, 0, 1), (-1, 0, 1), (0, 1, 1), (0, -1, 1)]
        elif crystal_type == "bcc":
            # BCC: {110}, {112}, and {123} planes with <111> directions
            planes = [(1, 1, 0), (1, 1, 2), (1, 2, 3)]
            directions = [(1, 1, 1), (-1, 1, 1), (1, -1, 1), (-1, -1, 1)]
        elif crystal_type == "sc":
            # SC: {100} planes with <100> directions
            planes = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
            directions = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        else:
            raise ValueError("Unsupported cubic crystal type: " + crystal_type)

    elif crystal_system == 'hexagonal' or crystal_system == 'trigonal':
        # Basal planes and <a> directions
        planes = [(0, 0, 1)]
        directions = [(1, 0, 0), (0, 1, 0), (-1, -1, 0)]

    elif crystal_system == 'tetragonal':
        # {001} planes with <100> directions
        planes = [(0, 0, 1)]
        directions = [(1, 0, 0), (0, 1, 0)]

    elif crystal_system == 'orthorhombic':
        planes = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        directions = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

    elif crystal_system == 'monoclinic':
        planes = [(0, 1, 0)]
        directions = [(1, 0, 0), (0, 0, 1)]

    elif crystal_system == 'triclinic':
        # For triclinic crystals, slip systems are not standardized
        raise NotImplementedError("Slip system generation for triclinic systems is not implemented.")

    else:
        raise ValueError("Unsupported crystal system.")

    # Generate slip systems
    slip_systems_list = []
    for plane in planes:
        for direction in directions:
            plane_vector = np.array(plane)
            direction_vector = np.array(direction)
            # Ensure the direction is within the plane
            if np.dot(plane_vector, direction_vector) == 0:
                slip_systems_list.append((plane, direction))

    # Remove duplicates
    slip_systems = list(set(slip_systems_list))
    return slip_systems, crystal_system



def rotate(pos, alpha, beta, gamma):
    """
    Rotates the lattice vectors in a POSCAR file by the specified Euler angles.

    Args:
        pos (str): Path to the POSCAR file.
        alpha (float): Rotation angle around the Z-axis in degrees.
        beta (float): Rotation angle around the Y-axis in degrees.
        gamma (float): Rotation angle around the X-axis in degrees.
    """
    import numpy as np
    import math

    # Read the POSCAR file
    with open(pos, 'r') as f:
        lines = f.readlines()

    # Extract lattice vectors from lines 3 to 5
    lattice_vectors = np.array([list(map(float, line.strip().split())) for line in lines[2:5]])

    # Ensure lattice vectors are 3D
    if lattice_vectors.shape != (3, 3):
        raise ValueError("Lattice vectors should be a 3x3 matrix.")

    # Convert angles from degrees to radians
    alpha_rad, beta_rad, gamma_rad = np.radians([alpha, beta, gamma])

    # Rotation matrices
    Rz = np.array([[math.cos(alpha_rad), -math.sin(alpha_rad), 0],
                   [math.sin(alpha_rad),  math.cos(alpha_rad),  0],
                   [0,                    0,                   1]])
    Ry = np.array([[ math.cos(beta_rad), 0, math.sin(beta_rad)],
                   [0,                   1, 0],
                   [-math.sin(beta_rad), 0, math.cos(beta_rad)]])
    Rx = np.array([[1, 0,                   0],
                   [0, math.cos(gamma_rad), -math.sin(gamma_rad)],
                   [0, math.sin(gamma_rad),  math.cos(gamma_rad)]])

    # Combined rotation matrix
    R = Rz @ Ry @ Rx  # Note: The order of multiplication matters

    # Apply rotation to lattice vectors
    rotated_vectors = lattice_vectors @ R.T

    # Write the rotated lattice vectors back to a new POSCAR file
    with open("POSCAR_rotated", 'w') as f:
        # Write the first two lines unchanged
        f.writelines(lines[:2])
        # Write the rotated lattice vectors with high precision
        for vec in rotated_vectors:
            f.write('  {:.16f}  {:.16f}  {:.16f}\n'.format(*vec))
        # Write the rest of the file unchanged
        f.writelines(lines[5:])

    print("Structure rotated and saved as 'POSCAR_rotated'.")


def calculate_schmid_factor(normal_vector, slip_direction, stress_direction):
    """
    Calculate the Schmid factor for a given slip system.
    :param normal_vector: Normal vector to the slip plane (as a list or tuple of 3 elements)
    :param slip_direction: Slip direction (as a list or tuple of 3 elements)
    :param stress_direction: Direction of the applied stress (as a list or tuple of 3 elements)
    :return: Schmid factor for the given slip system
    """

    # Normalize the vectors
    normal_vector = np.array(normal_vector)
    normal_vector = normal_vector / np.linalg.norm(normal_vector)

    slip_direction = np.array(slip_direction)
    slip_direction = slip_direction / np.linalg.norm(slip_direction)

    stress_direction = np.array(stress_direction)
    stress_direction = stress_direction / np.linalg.norm(stress_direction)

    # Calculate the Schmid factor
    m = np.dot(normal_vector, stress_direction) * np.dot(slip_direction, stress_direction)

    return abs(m)  # Schmid factor is usually taken as absolute value


