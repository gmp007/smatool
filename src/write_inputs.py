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

import os
from datetime import datetime

def write_default_input(cwd,dim,code_type):

    if dim =="1D":
        kpoints_content = """# Step: Static Calculation
0
G
4 4 2
0 0 0
# Step: Dynamical Calculation
0
G
4 4 2
0 0 0
"""
    elif dim =="2D":
        kpoints_content = """# Step: Static Calculation
0
G
4 4 1
0 0 0
# Step: Dynamical Calculation
0
G
4 4 1
0 0 0
"""
    elif dim =="3D":
        kpoints_content = """# Step: Static Calculation
0
G
4 4 4
0 0 0
# Step: Dynamical Calculation
0
G
4 4 4
0 0 0
"""
    # Write KPOINTS content to a file
    with open(os.path.join(cwd, "KPOINTS-sd"), "w") as kpoints_file:
        kpoints_file.write(kpoints_content)
    if code_type == "VASP":
    # INCAR content
        if dim in ["1D","3D"]:
            incar_content = """# Step: DFT Optimization
PREC    = Accurate
ENCUT   = 500
EDIFF   = 1e-6
EDIFFG  = -0.001
IBRION  = 2
ISIF    = 3
ISYM    = 2
NSW     = 300
ISMEAR  = 0
SIGMA   = 0.1
POTIM   = 0.1
PSTRESS = 0.001
LREAL   = False
NPAR    = 4
NSIM    = 4
ALGO    = Normal
IALGO   = 48
ISTART  = 0
LCHARG  = .FALSE.
LPLANE  = .TRUE.
LWAVE   = .FALSE.
# Step: DFT Yield Strength
PREC    = Accurate
ENCUT   = 600
EDIFF   = 1e-6
EDIFFG  = 1e-5
IBRION  = 2
ISIF    = 2
ISYM    = 2
ISMEAR  = 0
SIGMA   = 0.1
LREAL	= False
NPAR    = 11
NSW     = 300
ISTART  = 0
LCHARG  = .FALSE.
LPLANE  = .FALSE.
LWAVE   = .FALSE.
# Step: MD Nostrain
ENCUT   = 600
EDIFF   = 1E-6
EDIFFG  = 1e-5
ALGO    = Normal
IALGO   = 48
MAXMIX  = 40
IBRION  = 0
NBLOCK  = 1
KBLOCK  = 10
POTIM   = 1
ISYM    = 0
ISIF    = 3
TEBEG   = 300
LREAL   = False
NSW     = 100
PREC    = Normal
ISTART  = 0
ISMEAR  = 2
SIGMA   = 0.2
NPAR    = 11
NCORE   = 1
NSIM    = 4
NELMIN  = 4
NWRITE  = 0
LCHARG  = .FALSE.
LPLANE  = .TRUE.
LWAVE   = .FALSE.
IWAVPR  = 11
# Step: MD Yield Strength
ENCUT   = 600
EDIFF   = 1E-6
EDIFFG  = 1e-5
ALGO    = Normal
IALGO   = 48
MAXMIX  = 40
IBRION  = 0
NBLOCK  = 1
KBLOCK  = 10
POTIM   = 1
ISYM    = 0
ISIF    = 2
SMASS   = 2
MDALGO  = 2
NSW     = 500
TEBEG   = 300
PSTRESS = 0.0001
LREAL   = False
PREC    = Normal
ISTART  = 0
ISMEAR  = 2
SIGMA   = 0.2
NPAR    = 11
NCORE   = 1
NSIM    = 4
NELMIN  = 4
NWRITE  = 0
LCHARG  = .FALSE.
LPLANE  = .TRUE.
LWAVE   = .FALSE.
IWAVPR  = 11
"""
        elif dim == "2D":
            incar_content = """# Step: DFT Optimization
PREC    = Accurate
ENCUT   = 500
EDIFF   = 1e-6
EDIFFG  = -0.001
IBRION  = 2
ISIF    = 4
ISYM    = 2
NSW     = 300
ISMEAR  = 0
SIGMA   = 0.1
POTIM   = 0.1
PSTRESS = 0.001
LREAL   = False
NPAR    = 11
NSIM    = 4
ALGO    = Normal
IALGO   = 48
ISTART  = 0
LVDW	= True
IVDW	= 12
LCHARG  = .FALSE.
LPLANE  = .TRUE.
LWAVE   = .FALSE.
# Step: DFT Yield Strength
PREC    = Accurate
ENCUT   = 600
EDIFF   = 1e-6
EDIFFG  = 1e-5
IBRION  = 2
ISIF    = 2
ISYM    = 2
ISMEAR  = 0
SIGMA   = 0.1
LREAL	= False
NPAR    = 4
NSW     = 300
ISTART  = 0
LVDW    = True
IVDW    = 12
LCHARG  = .FALSE.
LPLANE  = .FALSE.
LWAVE   = .FALSE.
# Step: MD Nostrain
ENCUT   = 600
EDIFF   = 1E-6
EDIFFG  = 1e-5
ALGO    = Normal
IALGO   = 48
MAXMIX  = 40
IBRION  = 0
NBLOCK  = 1
KBLOCK  = 10
POTIM   = 1
ISYM    = 0
ISIF    = 4
TEBEG   = 300
LREAL   = False
NSW     = 100
PREC    = Normal
ISTART  = 0
ISMEAR  = 2
SIGMA   = 0.2
NPAR    = 11
NCORE   = 1
NSIM    = 4
NELMIN  = 4
NWRITE  = 0
LCHARG  = .FALSE.
LPLANE  = .TRUE.
LWAVE   = .FALSE.
IWAVPR  = 11
LVDW    = True
IVDW    = 12
# Step: MD Yield Strength
ENCUT   = 600
EDIFF   = 1E-6
EDIFFG  = 1e-5
ALGO    = Normal
IALGO   = 48
MAXMIX  = 40
IBRION  = 0
NBLOCK  = 1
KBLOCK  = 10
POTIM   = 1
ISYM    = 0
ISIF    = 2
SMASS   = 2
MDALGO  = 2
NSW     = 200
TEBEG   = 300
PSTRESS = 0.0001
LREAL   = False
PREC    = Normal
ISTART  = 0
ISMEAR  = 2
SIGMA   = 0.2
NPAR    = 11
NCORE   = 1
NSIM    = 4
NELMIN  = 4
NWRITE  = 0
LCHARG  = .FALSE.
LPLANE  = .TRUE.
LWAVE   = .FALSE.
IWAVPR  = 11
LVDW    = True
IVDW    = 12
"""
    elif code_type == "QE":
        if dim in ["1D","3D"]:
            qe_content = """{
  "steps": [
    {
      "name": "DFT Optimization",
      "data": {
        "input_data": {
          "control": {
            "calculation": "vc-relax",
            "restart_mode": "from_scratch",
            "pseudo_dir": "base_path",
            "tstress": true,
            "tprnfor": true,
            "forc_conv_thr": 0.00001,
            "outdir": "./OPT"
          },
          "system": {
            "ecutwfc": 60,
            "ecutrho": 600,
            "occupations": "smearing",
            "smearing": "mp",
            "degauss": 0.01
          },
          "electrons": {
            "conv_thr": 1e-8
          },
          "cell": {
            "cell_dofree": "all",
            "press": 0.0,
            "press_conv_thr": 0.5
          },
          "ions": {
            "tolp": 1.0e-4  
        },
        "pseudopotentials": "pseudopotentials",
        "kpts": "kpts"

        }
      }
    },
    {
      "name": "DFT Yield Strength",
      "data": {
        "input_data": {
          "control": {
            "calculation": "relax",
            "restart_mode": "from_scratch",
            "pseudo_dir": "base_path",
            "tstress": true,
            "tprnfor": true,
            "forc_conv_thr": 0.00001,
            "outdir": "./OPT"
          },
          "system": {
            "ecutwfc": 60,
            "ecutrho": 600,
            "occupations": "smearing",
            "smearing": "mp",
            "degauss": 0.01
          },
          "electrons": {
            "conv_thr": 1e-8
          },
          "cell": {
            "cell_dofree": "all",
            "press": 0.0,
            "press_conv_thr": 0.5
          },
          "ions": {
            "tolp": 1.0e-4
        },
        "pseudopotentials": "pseudopotentials",
        "kpts": "kpts"

        }
      }
    },
    {
      "name": "MD Nostrain",
      "data": {
        "input_data": {
          "control": {
            "calculation": "vc-md",
            "restart_mode": "from_scratch",
            "pseudo_dir": "base_path",
            "tstress": true,
            "tprnfor": true,
            "forc_conv_thr": 0.00001,
            "outdir": "./OPT"
          },
          "system": {
            "ecutwfc": 60,
            "ecutrho": 600,
            "occupations": "smearing",
            "smearing": "mp",
            "degauss": 0.01
          },
          "electrons": {
            "conv_thr": 1e-8
          },
          "cell": {
            "cell_dofree": "all",
            "press": 0.0,
            "press_conv_thr": 0.5
          },
          "ions": {
            "tolp": 1.0e-4,
            "tempw": 300
        },
        "pseudopotentials": "pseudopotentials",
        "kpts": "kpts"

        }
      }
    },
    {
      "name": "MD Yield Strength",
      "data": {
        "input_data": {
          "control": {
            "calculation": "md",
            "restart_mode": "from_scratch",
            "pseudo_dir": "base_path",
            "tstress": true,
            "tprnfor": true,
            "forc_conv_thr": 0.00001,
            "outdir": "./OPT"
          },
          "system": {
            "ecutwfc": 60,
            "ecutrho": 600,
            "occupations": "smearing",
            "smearing": "mp",
            "degauss": 0.01
          },
          "electrons": {
            "conv_thr": 1e-8
          },
          "cell": {
            "cell_dofree": "all",
            "press": 0.0,
            "press_conv_thr": 0.5
          },
          "ions": {
            "tolp": 1.0e-4,
            "tempw": 300
        },
        "pseudopotentials": "pseudopotentials",
        "kpts": "kpts"

        }
      }
    }
  ]
}

"""
        else:
            qe_content = """{
  "steps": [
    {
      "name": "DFT Optimization",
      "data": {
        "input_data": {
          "control": {
            "calculation": "vc-relax",
            "restart_mode": "from_scratch",
            "pseudo_dir": "base_path",
            "tstress": true,
            "tprnfor": true,
            "forc_conv_thr": 0.00001,
            "outdir": "./OPT"
          },
          "system": {
            "ecutwfc": 60,
            "ecutrho": 600,
            "occupations": "smearing",
            "smearing": "mp",
            "degauss": 0.01,
            "vdw_corr": "dft-d"
          },
          "electrons": {
            "conv_thr": 1e-8
          },
          "cell": {
            "cell_dofree": "2Dshape",
            "press": 0.0,
            "press_conv_thr": 0.5
          },
          "ions": {
            "tolp": 1.0e-4
        },
        "pseudopotentials": "pseudopotentials",
        "kpts": "kpts"

        }
      }
    },
    {
      "name": "DFT Yield Strength",
      "data": {
        "input_data": {
          "control": {
            "calculation": "relax",
            "restart_mode": "from_scratch",
            "pseudo_dir": "base_path",
            "tstress": true,
            "tprnfor": true,
            "forc_conv_thr": 0.00001,
            "outdir": "./OPT"
          },
          "system": {
            "ecutwfc": 60,
            "ecutrho": 600,
            "occupations": "smearing",
            "smearing": "mp",
            "degauss": 0.01,
            "vdw_corr": "dft-d"
          },
          "electrons": {
            "conv_thr": 1e-8
          },
          "cell": {
            "cell_dofree": "2Dshape",
            "press": 0.0,
            "press_conv_thr": 0.5
          },
          "ions": {
            "tolp": 1.0e-4
        },
        "pseudopotentials": "pseudopotentials",
        "kpts": "kpts"

        }
      }
    },
    {
      "name": "MD Nostrain",
      "data": {
        "input_data": {
          "control": {
            "calculation": "vc-md",
            "restart_mode": "from_scratch",
            "pseudo_dir": "base_path",
            "tstress": true,
            "tprnfor": true,
            "forc_conv_thr": 0.00001,
            "outdir": "./OPT"
          },
          "system": {
            "ecutwfc": 60,
            "ecutrho": 600,
            "occupations": "smearing",
            "smearing": "mp",
            "degauss": 0.01,
            "vdw_corr": "dft-d"
          },
          "electrons": {
            "conv_thr": 1e-8
          },
          "cell": {
            "cell_dofree": "2Dshape",
            "press": 0.0,
            "press_conv_thr": 0.5
          },
          "ions": {
            "tolp": 1.0e-4,
            "tempw": 300
        },
        "pseudopotentials": "pseudopotentials",
        "kpts": "kpts"

        }
      }
    },
    {
      "name": "MD Yield Strength",
      "data": {
        "input_data": {
          "control": {
            "calculation": "md",
            "restart_mode": "from_scratch",
            "pseudo_dir": "base_path",
            "tstress": true,
            "tprnfor": true,
            "forc_conv_thr": 0.00001,
            "outdir": "./OPT"
          },
          "system": {
            "ecutwfc": 60,
            "ecutrho": 600,
            "occupations": "smearing",
            "smearing": "mp",
            "degauss": 0.01,
            "vdw_corr": "dft-d"
          },
          "electrons": {
            "conv_thr": 1e-8
          },
          "cell": {
            "cell_dofree": "2Dshape",
            "press": 0.0,
            "press_conv_thr": 0.5
          },
          "ions": {
            "tolp": 1.0e-4,
            "tempw": 300
        },
        "pseudopotentials": "pseudopotentials",
        "kpts": "kpts"

        }
      }
    }
  ]
}
"""


    # Write INCAR content to a file
    if code_type == "VASP":
        with open(os.path.join(cwd, "INCARs"), "w") as incar_file:
            incar_file.write(incar_content)
    elif code_type == "QE":
        with open(os.path.join(cwd, "qe_input.in"), "w") as qe_file:
            qe_file.write(qe_content)



def write_default_ystool_in(cwd):
    if not os.path.exists(os.path.join(cwd, "smatool.in")):
        
        ystool_in = """########################################
###  SMATool package input control   ###
########################################
#choose stress calculator: VASP/QE currently supported
code_type = vasp

#choose method: DFT(static) and MD(dynamic)
mode = dft

#structure file name with .cif or .vasp
structure_file = file.cif

#choose dimension: 1D/2D/3D
dimensional = 2D

# strain: start, end, and interval
strains = 0.01 0.8 0.01

#Choose b/w tensile (off), compressional (compress) or both (on) 
plusminus = off

#yieldpoint offset, e.g. 0.2%
yieldpoint_offset = 0.002

#stress components: Tensile_x Tensile_y Tensile_z Tensile_biaxial indent_strength ideal_strength ideal_strength Shear (xz for 1D and 3D and xy for 2D)
# you can writ more than one case by just listing them (no comma, just space)  
components = Tensile_x

#supercell for MD simulation
md_supercell = 1, 1, 1

#save structure files at each strain
save_struct = on

#molecular dynamics time step
md_timestep = 500

#slipping on; must be on/yes/1 to model slip system along strain (plane) and slip directions
slipon = off

#define strain (plane) and slip directions, uses Miller indices, e.g., 11, 0-1 for 2D
slip_parameters = 100, 111

#define indentation radius in unit of your cell. Only needed when performing indentation strength
indent_radius = 2.0

#explicit potential directory
potential_dir = /potential/

#Use saved data. postprocessing only
use_saved_data = false

#job submission command
job_submit_command = vasp/pw.x > log
"""
        # Write to "elastool.in" file
        with open(os.path.join(cwd, "smatool.in"), "w") as elin:
            elin.write(ystool_in)





def print_default_input_message_0():
    print("╔════════════════════════════════════════════════════════════════════════════════╗")
    print("║                                                                                ║")
    print("║                       ♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥                       ║")
    print("║                  ♥♥♥  Default smatool.in template generated.  ♥♥♥              ║")
    print("║                 ♥♥ Modify and rerun smatool -0 to generate other   ♥♥          ║")
    print("║                 ♥♥    important input files. Happy running :)    ♥♥            ║")
    print("║                       ♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥                       ║")
    print("║                                   Exiting...                                   ║")
    print("║                                                                                ║")
    print("╚════════════════════════════════════════════════════════════════════════════════╝")





def print_default_input_message_1():
    print("╔════════════════════════════════════════════════════════════════════════════════╗")
    print("║                                                                                ║")
    print("║                               ♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥                                  ║")
    print("║                           ♥♥                    ♥♥                             ║")
    print("║                        ♥ All default inputs written to files. ♥                ║")
    print("║                        ♥ Modify according to dimensionality ♥                  ║")
    print("║                        ♥ and other criteria ♥                                  ║")
    print("║                        ♥       Happy running :)        ♥                       ║")
    print("║                           ♥♥                    ♥♥                             ║")
    print("║                               ♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥                                  ║")
    print("║                                Exiting...                                      ║")
    print("║                                                                                ║")
    print("╚════════════════════════════════════════════════════════════════════════════════╝")


max_width = len("|WARNING: This is an empirical approx; validity needs to be checked !! |")

def print_line(ec_file,content, padding=1, border_char="|", filler_char=" "):
    content_width = int(max_width) - (2 * int(padding)) - 2  # Subtract 2 for the border characters
    content = content[:content_width]  # Ensure content doesn't exceed the width
    line = border_char + filler_char*padding + content.ljust(content_width) + filler_char*padding + border_char
    #print(line)  # Print the line to the console
    if ec_file:
        ec_file.write(line + "\n")
    else:
        print(line)    

        
        

def print_banner2(yield_file, elapsed_time,ec_file=None, filename_struct=None):
    current_time = datetime.now().strftime('%H:%M:%S')
    current_date = datetime.now().strftime('%Y-%m-%d')
    
    # Adjust result message based on the presence of filename_struct
    if filename_struct:
        result_msg = f"Results are written in {filename_struct}\nand structure information\nat each stress and strain saved\nin {yield_file}."
    else:
        result_msg = f"Structure information\nat each stress and strain saved\nin {yield_file}."

    finish_msg = f"Job finished at {current_time} on {current_date} using {elapsed_time:.2f} s"

    message = f"{result_msg}\n{finish_msg}"

    print_line(ec_file,'❤' * (max_width - 2), padding=0, border_char='❤', filler_char='❤')
    for line in message.split('\n'):
        centered_line = line.center(max_width - 4)
        print_line(ec_file,centered_line, padding=1, border_char='❤')
    print_line(ec_file,'❤' * (max_width - 2), padding=0, border_char='❤', filler_char='❤')


def print_banner(version,component,code_type,method,ec_file=None):
    current_time = datetime.now().strftime('%H:%M:%S')
    current_date = datetime.now().strftime('%Y-%m-%d')
    conclusion_msg = f"Calculations started at {current_time} on {current_date}"
    component = component.upper() 

    message = f"Results for {component}\nusing\nSMATool Version: {version}\n {code_type} code is used as a calculator\nto perform {method} simulations\n{conclusion_msg}"

    max_width = 80  # Define the maximum width for the banner

    print_line(ec_file,'❤' * (max_width - 2), padding=0, border_char='❤', filler_char='❤')
    for line in message.split('\n'):
        centered_line = line.center(max_width - 4)
        print_line(ec_file,centered_line, padding=1, border_char='❤')
    print_line(ec_file,'❤' * (max_width - 2), padding=0, border_char='❤', filler_char='❤')




def print_boxed_message(ec_file=None):
    header_footer = "+" + "-" * 78 + "+"
    spacer = "| " + " " * 76 + " |"

    # List of lines to be printed
    lines = [
        (" * CITATIONS *", True),
        ("If you have used SMATool in your research, PLEASE cite:", False),
        ("", False),  # Space after the above line
        ("SMATool: Strength of materials analysis toolkit, ", False),
        ("C.E. Ekuma, ", False),
        ("Computer Physics Communications 300, 109189, (2024)", False),
        ("", False),

        ("", False),  # Blank line for separation
        ("Material strength analyzer: Python-based computational toolkit", False),
        ("for tensile and indent strength of materials, C.E. Ekuma,", False),
        ("www.github.com/gmp007/smatool", False)
    ]

    def output_line(line):
        if ec_file:
            ec_file.write(line + "\n")
        else:
            print(line)

    output_line(header_footer)
    
    for line, underline in lines:
        centered_line = line.center(76)
        output_line("| " + centered_line + " |")
        
        if underline:
            underline_str = "-" * len(centered_line)
            output_line("| " + underline_str.center(76) + " |")

    # Print footer of the box
    output_line(header_footer)
    

