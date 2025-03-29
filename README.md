# SMATool - Strength of Materials Analysis Toolkit 


**`SMATool`** — Strength of Materials Analysis Toolkit is a comprehensive computational package developed to evaluate the tensile, shear, biaxial, and indentation strength (Vickers' hardness) of materials at zero temperature (using density functional theory) and finite temperature (using ab-initio molecular dynamics). Integrating advanced algorithms and leveraging the capabilities of well-established electronic structure codes like **VASP** and **Quantum ESPRESSO** as the stress calculator, SMATool provides a robust platform for analyzing the strength of 1D, 2D, and 3D materials.

---

## ⚙️ Mechanical Output Summary

Below is a sample output summary for bulk cubic BN obtained using the VASP electronic structure code as a calculator. 


| Property                                               | Value              | Units       |
|--------------------------------------------------------|--------------------|-------------|
| Ultimate Strain (ε)                                    | 0.470              | —           |
| Ultimate Strength (σ<sub>s</sub>)                      | 190.350            | GPa         |
| Yield Strength (σ<sub>y</sub>)                         | 63.635             | GPa         |
| Stiffness (E)                                           | 719.792            | GPa         |
| Stiffness (K) (E-V)                                     | 934.371            | GPa         |
| Stiffness (K) (P-V)                                     | 735.995            | GPa         |
| Bulk Modulus Pressure Derivative (B<sub>p</sub>)        | 3.028              | —           |
| Thickness of Material                                   | —                  | —           |
| Resilience (U<sub>r</sub>)                             | 2.864              | J/m³        |
| Toughness (U<sub>t</sub>)                              | 65.507             | J/m³        |
| Energy Storage Capacity (E<sub>c</sub>)                | 5.592              | MJ/L        |
| Energy Storage Capacity (E<sub>c</sub>)                | 449.214            | Wh/kg       |
| Ductility (%D)                                          | 47.000             | %           |
| Ductility Index (%DI)                                  | 422.222            | %           |
| Poisson's Ratio (ν)                                     | 0.337              | —           |
| Shear Modulus (G)                                       | 269.181            | GPa         |
| Longitudinal Sound Velocity (V<sub>l</sub>)            | 14.591             | km/s        |
| Shear Sound Velocity (V<sub>s</sub>)                   | 8.824              | km/s        |
| Mean Sound Velocity (V<sub>m</sub>)                    | 9.904              | km/s        |
| P-wave Velocity (V<sub>p</sub>)                        | 17.796             | km/s        |
| Debye Temperature (T<sub>D</sub>)                      | 1625.905           | K           |
| First Lamé Parameter (λ)                               | 556.541            | GPa         |
| Second Lamé Parameter (μ)                              | 269.181            | GPa         |
| Grüneisen Parameter (γ)                                 | 2.028              | —           |
| Min Thermal Conductivity (Slack)                        | 0.000              | W/m·K       |
| Min Thermal Conductivity (Clarke)                       | 1.318              | W/m·K       |
| Min Thermal Conductivity (Cahill)                       | 6.976              | W/m·K       |
| Integrated Min Thermal Conductivity (Slack)             | 2.674              | W/m·K       |


## Importance of SMATool
The significance of SMATool lies in its ability to provide accurate and efficient strength parameters, which are crucial in material design and analysis. Understanding these mechanical properties helps in predicting material behavior under different loading conditions and in different environments, contributing to the development of materials with optimized performance characteristics.

## 🔍 Why SMATool?

The **importance of SMATool** lies in its ability to provide *accurate*, *automated*, and *efficient* evaluation of strength parameters, which are critical for material design and failure prediction. These insights enable:

- **Predictive materials design**
- **Optimization of mechanical performance**
- **Failure analysis under complex loading conditions**
- **Cross-validation with experiments and machine learning models**


## Types of Strengths within SMATool
1. *Uniaxial Tensile Strength:* The maximum stress a material can withstand while being stretched or pulled before failure. Tensile strength is crucial in applications where materials are subject to tension. The tensile strength can be computed along the various directions: x, y, and z, and along the critical plane directions. 

2. *Biaxial Strength:* The strength of a material when subjected to biaxial stress conditions. This is important for materials used in applications where they are exposed to stresses in two perpendicular directions.

3. *Shear Strength:* This measures a material's ability to resist shear loads, which is vital in understanding how materials will behave when forces are applied parallel to a material's surface. The shear strength can be determined along the critical plane directions and slip directions. 


4. *Indentation Strength:* The measure of the resistance of a material to deformation under localized compressive forces. This can be used to quantify the Vickers' hardness of a material.

Additionally, the SMATool package computes the yield strength for each of the strengths above, and the corresponding energy storage capacity (in MJ/L and Wh/Kg) of the material at the ultimate strain. 

## 🧮 SMATool Stress Calculators

SMATool is powered by:
- **VASP (Vienna Ab initio Simulation Package)**
- **Quantum ESPRESSO**
- Integrated via the **ASE (Atomic Simulation Environment)**

These calculators ensure:
- High-fidelity stress-strain data
- Automation of stress tensor tracking
- Consistency across different material classes

💡 *We welcome collaborative implementation of other codes (e.g., ABINIT, CASTEP, SIESTA, etc.) upon request.*


## Installation
**SMATool** offers straightforward installation options suitable for various user preferences as explained below. We note that in all the installation options, all the libraries and dependables are automatically determined and installed alongside the SMATool. 

1. **Using pip**: Our recommended way to install the **SMATool** package is using pip. 
   - Quickly install the latest version of the SMATool package with pip by executing: 
     ```
     pip install -U smatool
     ```

2. **From Source Code**:
   - Alternatively, users can download the source code with:
     ```
     git clone [git@github.com:gmp007/smatool.git]
     ```
   - Then, install SMATool by navigating to the master directory and running:
     ```
     pip install .
     ```

3. **Installation via setup.py**:
   - SMATool can also be installed using the `setup.py` script:
     ```
     python setup.py install [--prefix=/path/to/install/]
     ```
   - The optional `--prefix` argument is useful for installations in environments like shared High-Performance Computing (HPC) systems, where administrative privileges might be restricted.
   - Please note that while this method remains supported, its usage is gradually declining in favor of more modern installation practices. We only recommend this installation option where standard installation methods like `pip` are not applicable.
   
## Usage and Running SMATool

The best way to learn how to use the SMATool package is to start with the provided examples folder. The key steps for initializing SMATool follows:

1. **Create a Calculation Directory**:
   - Start by creating a directory for your calculations.
   - Run `smatool -0` to generate main input template of the SMATool, which is the `smatool.in`.

2. **Modify Input Files**:
   - Customize the generated files according to your project's requirements, choose the code type between VASP and QE, and specify the directory of your potential files. 

3. **Initialize the Job**:
   - Execute `smatool` to begin the calculation process.

4. **Understanding SMATool Options**:
   - The main input file `smatool.in` includes descriptive text for each flag, making it user-friendly.

## Citing SMATool
If you have used the SMATool package in your research, please cite:
  - [SMATool: An automated toolkit for strength of materials](https://doi.org/10.1016/j.cpc.2024.109189) - 

@article{Ekuma2024,
  title = {SMATool: Strength of Materials Analysis Toolkit},
  journal = {Computer Programs in Physics},
  volume = {300},
  pages = {109189},
  year = {2024},
  doi = {https://doi.org/10.1016/j.cpc.2024.109189},
  url = {https://doi.org/10.1016/j.cpc.2024.109189},
  author = {Chinedu E. Ekuma}
}

## SMATool Utility
VASP electronic structure calculations come with proprietary pseudopotentials included in the package. In contrast, Quantum Espresso (QE) is open source, offering a variety of sources for obtaining pseudopotentials. While no specific pseudopotential database is officially recommended for QE, we prefer the norm-conserving [Pdojo pseudopentials](http://www.pseudo-dojo.org/). The SMATool computational toolkit includes an automated utility package, `qepotential`, located in the utility folder. This tool automates the generation of pseudopotentials from the Pdojo website for all materials required to run the SMATool package with QE as the calculator. Users can also specify custom requirements in the `pseudo.info` input file. The SMATool utility package saves the potentials of various elements as `element_pdojo.upf` in the `qe_potentials` folder.


## 🚀 Coming Soon

- **Strain rate effects** for dynamic loading
- **Fracture mechanics** (including crack propagation modeling)
- Graphical tools for **stress-strain visualization**
- Plug-and-play **ML models** for property prediction


## Contact Information
Please if you find a bug, want to extend the SMATool computational toolkit, or have an idea that you would want us to incorporate, please reach out to us. Our team is dedicated to supporting your work and enhancing capabilities and efficiency of the SMATool package.

Feel free to contact us via email:
- [cekuma1@gmail.com](mailto:cekuma1@gmail.com)

Your feedback and questions are invaluable to us, and we look forward to hearing from you.

## License

This project is licensed under the GNU GPL version 3 - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- This work was supported by the U.S. Department of Energy, Office of Science, Basic Energy Sciences under Award DOE-SC0024099.

---
