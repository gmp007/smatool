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



class ShearTensileConstraint:
    """Custom constraint for shear and tensile stress."""
    def __init__(self, stress_component=None, dim="3D"):
        """
        stress_component: Direction for tensile and shear constraint.
        dim: Dimensional.
        """
        self.stress_component = stress_component
        self.dim =dim


    def adjust_positions(self, atoms, newpositions):
        # This method is required by ASE but does nothing in this case
        pass

    def index_shuffle(self, atoms, newindices):
        # This method is required by ASE but does nothing for this constraint
        pass
        
    def get_removed_dof(self, atoms):
        # This method can return the number of degrees of freedom removed by the constraint if applicable
        # For now, let's return 0 to indicate no explicit removal of degrees of freedom
        return 0
        
        
    def adjust_forces(self, atoms, forces):
        # Handle tensile constraints
        if self.stress_component == "Tensile_x":
            forces[:, 0] = 0  # Zero out x-component for tensile along xx
        elif self.stress_component == "Tensile_y":
            forces[:, 1] = 0  # Zero out y-component for tensile along yy
        elif self.stress_component == "Tensile_z":
            forces[:, 2] = 0  # Zero out z-component for tensile along zz
        elif self.stress_component in ["Shear", "indent_strength"]:
            if self.dim in ["3D","1D"]:
                forces[:, 1] = 0  # Zero out y-component for shear in xz plane
            else:
                forces[:, 2] = 0  # Zero out z-component for shear in xy plane
        elif self.stress_component == "Tensile_biaxial":
            forces[:, 0] = 0  # Zero out x-component for tensile along xx
            forces[:, 1] = 0  # Zero out y-component for tensile along yy
        else:
            print(f"Stress component '{self.stress_component}' not implemented.")            

