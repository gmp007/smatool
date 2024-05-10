# SMATool - Strength of Materials Analysis Toolkit 


**`SMATool`** -- Strength of Materials Analysis Toolkit is a comprehensive computational package developed to evaluate the tensile, shear, biaxial, and indentation strength (Vickers' hardness) of materials at zero temperature (using density functional theory) and finite temperature (using ab-initio molecular dynamics). Integrating advanced algorithms and leveraging the capabilities of well-established electronic structure codes like VASP and QE as the stress calculator, SMATool provides a robust platform for analyzing the strength of 1D, 2D, and 3D materials.

## Importance of SMATool
The significance of SMATool lies in its ability to provide accurate and efficient strength parameters, which are crucial in material design and analysis. Understanding these mechanical properties helps in predicting material behavior under different loading conditions and in different environments, contributing to the development of materials with optimized performance characteristics.

## Types of Strengths within SMATool
1. *Uniaxial Tensile Strength:* The maximum stress a material can withstand while being stretched or pulled before failure. Tensile strength is crucial in applications where materials are subject to tension. The tensile strength can be computed along the various directions: x, y, and z, and along the critical plane directions. 

2. *Biaxial Strength:* The strength of a material when subjected to biaxial stress conditions. This is important for materials used in applications where they are exposed to stresses in two perpendicular directions.

3. *Shear Strength:* This measures a material's ability to resist shear loads, which is vital in understanding how materials will behave when forces are applied parallel to a material's surface. The shear strength can be determined along the critical plane directions and slip directions. 


4. *Indentation Strength:* The measure of the resistance of a material to deformation under localized compressive forces. This can be used to quantify the Vickers' hardness of a material.

Additionally, the SMATool package computes the yield strength for each of the strengths above, and the corresponding energy storage capacity (in MJ/L and Wh/Kg) of the material at the ultimate strain. 

## SMATool Stress Calculators
SMATool utilizes VASP (Vienna Ab initio Simulation Package) and QE (Quantum Espresso) as its calculator, ensuring high accuracy and reliability. These powerful electronic structure codes, combined with the ASE (Atomic Simulation Environment) backend, enable detailed and precise material analysis. If you're interested in other electronic structure codes such as the stress calculation, we will be more than happy to collaborate in such implementations.


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
   - Run `smatool -0` to generate the main input template of the SMATool, which is the `smatool.in`.

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
  journal = {Computer Physics Communications},
  volume = {300},
  pages = {109189},
  year = {2024},
  doi = {https://doi.org/10.1016/j.cpc.2024.109189},
  url = {https://www.sciencedirect.com/science/article/abs/pii/S0010465524001127},
  author = {Chinedu Ekuma}
}

## SMATool Utility
ASP electronic structure calculations come with proprietary pseudopotentials included in the package. In contrast, Quantum Espresso (QE) is open source, offering a variety of sources for obtaining pseudopotentials. While no specific pseudopotential database is officially recommended for QE, we prefer the norm-conserving [Pdojo pseudopentials](http://www.pseudo-dojo.org/). The SMATool computational toolkit includes an automated utility package, `qepotential`, located in the utility folder. This tool automates the generation of pseudopotentials from the Pdojo website for all materials required to run the SMATool package with QE as the calculator. Users can also specify custom requirements in the `pseudo.info` input file. The SMATool utility package saves the potentials of various elements as `element_pdojo.upf` in the `qe_potentials` folder.

## Got Questions
To join the Google user group, post your questions, and see if your inquiries have already been addressed, you can visit the [SMA Tools User Group](https://groups.google.com/g/smatools/) on Google Groups. This platform allows for interactive discussions and access to previously answered questions, facilitating a comprehensive support community.


## Contact Information
Please if you find a bug, want to extend the SMATool computational toolkit, or have an idea that you would want us to incorporate, please reach out to us. Our team is dedicated to supporting your work and enhancing the capabilities and efficiency of the SMATool package.

Feel free to contact us via email:
- [cekuma1@gmail.com](mailto:cekuma1@gmail.com)

Your feedback and questions are invaluable to us, and we look forward to hearing from you.

## License

This project is licensed under the GNU GPL version 3 - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- This work was supported by the U.S. Department of Energy, Office of Science, Basic Energy Sciences under Award DOE-SC0024099.

---
