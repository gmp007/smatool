import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
import sys
from write_inputs import write_default_input

class FileHandler:
    def __init__(self, tmp_path):
        self.tmp_path = tmp_path

    def save_path(self, pp_path_entry):
        base_path = self.tmp_path
        try:
            with open(os.path.join(base_path, "pp_path.txt"), "w") as file:
                file.write(pp_path_entry.get() + "\n")
        except FileNotFoundError:
            pass

    def load_config(self, pp_path_entry):
        base_path = self.tmp_path
        try:
            with open(os.path.join(base_path, "pp_path.txt"), "r") as file:
                lines = file.readlines()
                if len(lines) >= 1:
                    pp_path_entry.insert(0, lines[0].strip())
        except FileNotFoundError:
            pass
            
class UFCFileGenerator:
    def __init__(self, tab, file_handler):
        self.file_handler = file_handler
        self.create_canvas(tab)
        self.create_form()

    def load_config(self):
        self.file_handler.load_config(self.pp_path_entry)

    def create_canvas(self, tab):
        def configure_scroll_region(event):
            canvas.configure(scrollregion=canvas.bbox("all"))

        canvas = tk.Canvas(tab)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        frame = tk.Frame(canvas)
        canvas.create_window((0, 0), window=frame, anchor=tk.NW)

        v_scrollbar = tk.Scrollbar(tab, orient=tk.VERTICAL, command=canvas.yview, bg='black')
        v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.config(yscrollcommand=v_scrollbar.set)

        h_scrollbar = tk.Scrollbar(tab, orient=tk.HORIZONTAL, command=canvas.xview, bg='black')
        h_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)
        canvas.config(xscrollcommand=h_scrollbar.set)

        canvas.bind_all("<MouseWheel>", lambda event: canvas.yview_scroll(-1 * int(event.delta / 120), "units"))
        canvas.bind_all("<Shift-MouseWheel>", lambda event: canvas.xview_scroll(-1 * int(event.delta / 120), "units"))

        frame.bind("<Configure>", configure_scroll_region)
        self.frame = frame

    def create_form(self):
        self.create_button_frame()
        self.create_form_fields()
        self.create_additional_fields()
        self.load_config()  # Call load_config after pp_path_entry is defined

    def create_button_frame(self):
        button_frame = tk.Frame(self.frame)
        button_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=10)

        tk.Button(button_frame, bg='turquoise', text="Generate smatool.in", command=self.generate_smatoolin_file).grid(row=0, column=0, padx=5, pady=1, sticky="w")
        tk.Button(button_frame, bg='#ffff85', text="View/Edit smatool.in", command=self.open_sample_file).grid(row=0, column=1, padx=5, pady=1, sticky="w")
        tk.Button(button_frame, text="Generate KPOINTS and INPUT FILES", command=self.generate_kpoints_incar).grid(row=0, column=2, padx=5, pady=1, sticky="w")
        tk.Button(button_frame, text="help", command=lambda kw="smatool_file": display_help(kw, use_tk=True)).grid(row=0, column=3, padx=5, pady=1, sticky="w")
        tk.Button(button_frame, text="Exit", command=self.exit_application).grid(row=0, column=4, padx=5, pady=1, sticky="w")

    def create_form_fields(self):
        self.mode = tk.LabelFrame(self.frame, text="SMATool INPUT: ", font=("Helvetica", 14, "bold"))
        self.mode.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=8, pady=(3, 3))

        self.create_dropdown("code_type", "Choose code type:", ["vasp", "qe"], 0)
        self.create_dropdown("mode", "Choose mode:", ["dft", "md"], 1)
        self.create_file_entry("structure_file", "Structure file name:", 2, "Browse", self.browse_structure_file)
        self.create_dropdown("dimensional", "Choose dimension:", ["1D", "2D", "3D"], 3)
        #self.create_strain_entry("strains", "Strain: start, end, interval", 4)
        self.create_dropdown("plusminus", "Perform compression with tensile strain:", ["off", "on"], 5)
        #self.create_entry("yieldpoint_offset", "Yield point offset:", 6)
        
        self.create_multiselect("components", "Stress components:", ["Tensile_x", "Tensile_y", "Tensile_biaxial", "indent_strength", "Shear"], 7)
        #self.create_entry("md_supercell", "Supercell for MD simulation:", 8)
        self.create_dropdown("save_struct", "Save structure files:", ["on", "off"], 9)
        #self.create_entry("md_timestep", "MD time step:", 10)
        #self.create_entry("indent_radius", "Indentation radius:", 11, default_value="2.0")
        #self.create_entry("slip_parameters", "Slip parameters:", 12)
        #self.create_dropdown("slipon", "Enable slipping:", ["on", "off"], 13)
        self.create_file_entry("potential_dir", "Potential directory:", 14, "Browse", self.browse_potential_dir)
        self.create_dropdown("use_saved_data", "Use saved data:", ["false", "true"], 15)
        #self.create_entry("job_submit_command", "Job submission command:", 16)


    def create_parallel_command_field(self, row):
        self.command_label = tk.Label(self.mode, text="Job submission command:").grid(row=row, column=0, sticky='w', pady=1)
        self.command_entry = tk.Entry(self.mode, width=30, bg='white')
        self.command_entry.grid(row=row, column=1, columnspan=9, sticky='we', pady=1, padx=5)
        self.command_entry.insert(0, "mpirun -n 4 vasp_std > vasp.log")
        
        
    def create_strain_entry(self, var_name, label_text, row):
        label = tk.Label(self.mode, text=label_text).grid(row=row, column=0, sticky='w', pady=1)
        entry = tk.Entry(self.mode, width=55, bg='white')
        entry.grid(row=row, column=1, columnspan=9, sticky='w', pady=1)
        entry.insert(0, "0.01 0.8 0.01")
        setattr(self, var_name + "_entry", entry)

    def create_multiselect(self, var_name, label_text, options, row):
        label = tk.Label(self.mode, text=label_text).grid(row=row, column=0, sticky='w', pady=1)
        var = tk.StringVar()
        var.set(" ".join(options))  # Default to all options selected
        menu = tk.OptionMenu(self.mode, var, *options)
        menu.config(width=55)
        menu.grid(row=row, column=1, columnspan=9, sticky='w', pady=1)
        setattr(self, var_name + "_var", var)
        setattr(self, var_name + "_menu", menu)

    def create_dropdown(self, var_name, label_text, options, row, command=None):
        label = tk.Label(self.mode, text=label_text).grid(row=row, column=0, sticky='w', pady=1)
        var = tk.StringVar()
        var.set(options[0])
        menu = tk.OptionMenu(self.mode, var, *options, command=command)
        menu.config(width=10)
        menu.grid(row=row, column=1, columnspan=9, sticky='w', pady=1)
        setattr(self, var_name + "_var", var)
        setattr(self, var_name + "_menu", menu)

    def create_entry(self, var_name, label_text, row, default_value=None, width=55):
        label = tk.Label(self.mode, text=label_text).grid(row=row, column=0, sticky='w', pady=1)
        entry = tk.Entry(self.mode, width=width, bg='white')
        entry.grid(row=row, column=1, columnspan=9, sticky='w', pady=1)
        if default_value is not None:
            entry.insert(0, default_value)
        setattr(self, var_name + "_entry", entry)
        
        
    def create_file_entry(self, var_name, label_text, row, button_text, button_command):
        label = tk.Label(self.mode, text=label_text).grid(row=row, column=0, sticky='w', pady=1)
        entry = tk.Entry(self.mode, width=55, bg='white')
        entry.grid(row=row, column=1, columnspan=9, sticky='w', pady=1)
        button = tk.Button(self.mode, text=button_text, command=button_command)
        button.grid(row=row, column=10, padx=5, pady=10, sticky="e")
        setattr(self, var_name + "_entry", entry)
        setattr(self, var_name + "_button", button)

    def browse_structure_file(self):
        filepath = filedialog.askopenfilename(title="Select structure file (cif or vasp format)")
        if filepath:
            self.structure_file_entry.delete(0, tk.END)
            self.structure_file_entry.insert(0, filepath)

    def browse_potential_dir(self):
        filepath = filedialog.askdirectory(title="Select potential directory")
        if filepath:
            self.potential_dir_entry.delete(0, tk.END)
            self.potential_dir_entry.insert(0, filepath)

    def open_sample_file(self):
        InputFileViewer("smatool.in")

    def generate_smatoolin_file(self):
        filename = "smatool.in"
        self.list_keys = self.collect_form_data()
        with open(filename, 'w') as f:
            self.write_smatoolin_file(f)
        messagebox.showinfo("Success", f"{filename} has been generated.")
        self.file_handler.save_path(self.pp_path_entry)

    def collect_form_data(self):
        list_keys = {
            "code_type": self.code_type_var.get(),
            "mode": self.mode_var.get(),
            "structure_file": self.structure_file_entry.get(),
            "dimensional": self.dimensional_var.get(),
            "strains": self.strains_entry.get(),
            "plusminus": self.plusminus_var.get(),
            "yieldpoint_offset": self.yieldpoint_offset_entry.get(),
            "components": self.components_var.get(),
            "md_supercell": self.supercell_size_entry.get(),
            "save_struct": self.save_struct_var.get(),
            "md_timestep": self.md_timestep_entry.get(),
            "indent_radius": self.indent_radius_entry.get(), 
            "slip_parameters": ', '.join(entry.get() for entry in self.slip_parameters_entries),  # Updated this line
            "slipon": self.slipon_var.get(),
            "potential_dir": self.pp_path_entry.get(),
            "use_saved_data": self.use_saved_data_var.get(),
            "job_submit_command": self.command_entry.get()
        }
        return list_keys


    def write_smatoolin_file(self, file):
        print("###############################################################################", file=file)
        print("### The input file to control the calculation details of SMATool package     ###", file=file)
        print("###############################################################################\n", file=file)
        
        print("# Choose electronic structure code calculator: VASP and QE currently supported", file=file)
        print(f"code_type = {self.list_keys['code_type']}\n", file=file)
        
        print("# Choose method for stress-strain: DFT(static) and MD(dynamic)", file=file)
        print(f"mode = {self.list_keys['mode']}\n", file=file)
        
        print("# Structure file name with .cif or .vasp", file=file)
        print(f"struct_file = {self.list_keys['structure_file']}\n", file=file)
        
        print("# Choose the dimension of your material: 1D/2D/3D", file=file)
        print(f"dimensional = {self.list_keys['dimensional']}\n", file=file)
        
        print("# Strain: start, end, and interval", file=file)
        print(f"strains = {self.list_keys['strains']}\n", file=file)
        
        print("# Perform compression with tensile strain", file=file)
        print(f"plusminus = {self.list_keys['plusminus']}\n", file=file)
        
        print("# Yield point offset => For metals 0.2%", file=file)
        print(f"yieldpoint_offset = {self.list_keys['yieldpoint_offset']}\n", file=file)
        
        print("# Stress components: Tensile_x Tensile_y Tensile_z Tensile_biaxial ideal_strength Shear (xz for 1D and 3D and xy for 2D)", file=file)
        print("# You can write more than one case by listing them with space", file=file)
        print(f"components = {self.list_keys['components']}\n", file=file)
        
        print("# Supercell for MD simulation", file=file)
        print(f"md_supercell = {self.list_keys['md_supercell']}\n", file=file)
        
        print("# Save structure files at each strain", file=file)
        print(f"save_struct = {self.list_keys['save_struct']}\n", file=file)
        
        print("# Molecular dynamics time step", file=file)
        print(f"md_timestep = {self.list_keys['md_timestep']}\n", file=file)
        
        print("# define indentation radius in unit of your cell. Only needed when performing indentation strength", file=file)
        print(f"indent_radius = {self.list_keys['indent_radius']}\n", file=file)
        
        print("# Define strain direction, slip direction; for 2D, omit last Miller indices, e.g., 11, 0-1", file=file)
        print(f"slip_parameters = {self.list_keys['slip_parameters']}\n", file=file)
        
        print("# Slipping on; must be on/yes/1 once shear modeling is enabled", file=file)
        print(f"slipon = {self.list_keys['slipon']}\n", file=file)
        
        print("# Potential directory", file=file)
        print(f"potential_dir = {self.list_keys['potential_dir']}\n", file=file)
        
        print("# Use saved data; postprocessing to obtain yield strength", file=file)
        print(f"use_saved_data = {self.list_keys['use_saved_data']}\n", file=file)
        
        print("# Job submission command", file=file)
        print(f"job_submit_command = {self.list_keys['job_submit_command']}\n", file=file)


    def generate_kpoints_incar(self):
        cwd = os.getcwd()
        code_type = self.code_type_var.get().upper()
        dim = self.dimensional_var.get()
        print(f"Generating KPOINTS and INCAR files with code_type: {code_type}, dimension: {dim}, cwd: {cwd}")
        write_default_input(code_type, dim, cwd)
        messagebox.showinfo("Success", "KPOINTS and INCAR files have been generated.")


            
    def create_repeat_numbers(self, row):
        lbl = tk.Label(self.mode, text="Supercell for MD simulation:")
        lbl.grid(row=row, column=0, sticky='w', pady=1)
        self.supercell_size_entry = tk.Entry(self.mode, width=10, bg='white')
        self.supercell_size_entry.grid(row=row, column=1, columnspan=3, sticky='w', padx=5, pady=1)
        self.supercell_size_entry.insert(0, "1, 1, 1")


    def create_slip_parameters(self, row):
        lbl = tk.Label(self.mode, text="Slip parameters (e.g., 100, 111):")
        lbl.grid(row=row, column=0, sticky='w', pady=1)
        
        self.slip_parameters_entries = []
        default_values = ["100", "111"]
        for i, default_value in enumerate(default_values):
            entry = tk.Entry(self.mode, width=10, bg='white')
            entry.grid(row=row, column=i + 1, sticky='w', padx=5, pady=1)
            entry.insert(0, default_value)
            self.slip_parameters_entries.append(entry)



    def create_additional_fields(self):
        row = 18  # Adjust the starting row to avoid overlap with create_form_fields
        

        self.create_entry("strains", "Strain: start, end, interval", row, default_value="0.01 0.8 0.01",width=15)
        
        row += 1

        
        self.create_entry("yieldpoint_offset", "Yield point offset:", row, default_value="0.002",width=10)
        
        row += 1
        # Add supercell for MD simulation
        #self.create_entry("md_supercell", "Supercell for MD simulation:", row)
        self.create_repeat_numbers(row)
        
        row += 1
        # Add MD time step
        #self.create_entry("md_timestep", "MD time step:", row)
        self.create_entry("md_timestep", "MD time step:", row, default_value="500",width=10)
        
        #row += 1
        ## Add save structure files option
        #self.create_dropdown("save_struct", "Save structure files:", ["on", "off"], row)
        
        row += 1
        # Add Indentation radions
        self.create_entry("indent_radius", "Indentation radius:", row, default_value="2.0",width=10)
        #self.create_dropdown("indent_radius", "Choose indentation radius:", ["on", "off"], row)
        
        row += 1
        # Add slip parameters
        #self.create_entry("slip_parameters", "Slip parameters:", row)
        self.create_slip_parameters(row)
        
        row += 1
        # Add slip on option
        self.create_dropdown("slipon", "Enable slipping:", ["off", "on"], row)
        
        #row += 1
        ## Add potential directory field
        #self.create_file_entry("potential_dir", "Potential directory:", row, "Browse", self.browse_potential_dir)
        
        #row += 1
        ## Add use saved data option
        #self.create_dropdown("use_saved_data", "Use saved data:", ["true", "false"], row)
        
        #row += 1
        # Add job submission command field
        #self.create_entry("job_submit_command", "Job submission command:", row)
        
        row += 2
        # Add pp_path_entry field initialization here
        self.create_pp_path_field(row)


        row += 2
        self.create_parallel_command_field(row)
        
        
    def create_pp_path_field(self, row):
        path_pp_label = tk.Label(self.mode, text="Pseudopotential path:")
        path_pp_label.grid(row=row, column=0, sticky='w', pady=1)
        self.pp_path_entry = tk.Entry(self.mode, width=55, bg='white')
        self.pp_path_entry.grid(row=row, column=1, columnspan=9, sticky='w', pady=1)
        self.browse_pp_button = tk.Button(self.mode, text="Browse", command=self.browse_pp_path)
        self.browse_pp_button.grid(row=row, column=10, padx=5, pady=10, sticky="e")

    def browse_pp_path(self):
        filepath = filedialog.askdirectory(title="Pseudo Potential directory - location of POTCAR files for VASP")
        if filepath:
            self.pp_path_entry.delete(0, tk.END)
            self.pp_path_entry.insert(0, filepath)


    def exit_application(self):
        print("======================================================")
        print("= SMATool input files successfully written to file   =")
        print("=  Modify if necessary and proceed with runs Exiting =")
        print("=          Happy simulations Exiting ...             =")
        print("======================================================")
        self.frame.quit()
        self.frame.update()
        self.frame.destroy()
        


def display_help(keyword=None, use_tk=False):
    help_text = """---------------------------
SMATool Help Documentation
---------------------------

SMATool is a software package for stress and mechanical properties analysis. It supports VASP and QE for electronic structure calculations and allows for stress-strain analysis using DFT or MD methods.

Configuration and Input Parameters:
-----------------------------------

code_type: Choose electronic structure code calculator (vasp/qe).
mode: Choose method for stress-strain analysis (dft/md).
structure_file: Specify the structure file name with .cif or .vasp extension.
dimensional: Specify the dimension of your material (1D/2D/3D).
strains: Define the strain parameters (start, end, interval).
plusminus: Enable/disable compression with tensile strain (on/off).
yieldpoint_offset: Set the yield point offset (e.g., 0.002 for metals).
components: List the stress components to be analyzed (Tensile_x, Tensile_y, Tensile_biaxial, indent_strength, Shear).
md_supercell: Define the supercell dimensions for MD simulation, format: x, y, z (e.g., 1, 1, 1).
save_struct: Enable/disable saving structure files at each strain (on/off).
md_timestep: Set the molecular dynamics time step (default is 500).
indent_radius: Define the indentation radius in units of your cell (default is 2.0).
rotation: Enable/disable rotation along strain and slip directions (on/off).
slip_parameters: Define the strain direction and slip direction (e.g., 100, 111).
slipon: Enable/disable slipping (on/off).
potential_dir: Specify the potential directory for VASP or QE.
use_saved_data: Enable/disable using saved data for postprocessing (true/false).
job_submit_command: Define the job submission command.

Additional Information:
------------------------

- code_type: Choose between 'vasp' for VASP calculations or 'qe' for Quantum ESPRESSO calculations.
- mode: Choose 'dft' for Density Functional Theory calculations or 'md' for Molecular Dynamics simulations.
- strains: Define the strain range and intervals. For example, '0.01 0.8 0.01' means starting at 0.01, ending at 0.8 with an interval of 0.01.
- plusminus: Set 'on' to perform both tensile and compressive strains or 'off' for tensile only.
- components: Specify the stress components to analyze. Multiple components can be specified separated by space.
- md_supercell: Define the size of the supercell for MD simulations in the format x, y, z.
- save_struct: When set to 'on', structure files will be saved at each strain step.
- md_timestep: Set the timestep for MD simulations, typically in femtoseconds.
- indent_radius: Define the radius for indentation strength calculations. This is required when performing indentation tests.
- rotation: If set to 'on', the system will rotate along specified strain and slip directions.
- slip_parameters: Define the directions for strain and slip. For example, '100, 111'.
- slipon: Set to 'on' to enable slipping which is required for shear modeling.
- potential_dir: Specify the directory containing potential files for VASP or QE.
- use_saved_data: When set to 'true', the program will use previously saved data for postprocessing to obtain yield strength.
- job_submit_command: Provide the command used to submit jobs, such as 'mpirun -n 2 vasp_std > log'.

Example Configuration:
-----------------------

code_type = vasp
mode = dft
structure_file = cBN.cif
dimensional = 3D
strains = 0.01 0.8 0.01
plusminus = off
yieldpoint_offset = 0.002
components = Tensile_x Tensile_y Tensile_biaxial indent_strength Shear
md_supercell = 1, 1, 1
save_struct = on
md_timestep = 500
indent_radius = 2.0
rotation = off
slip_parameters = 100, 111
slipon = off
potential_dir = /home/user/codes/vasp/PBE
use_saved_data = true
job_submit_command = mpirun -n 2 vasp_std > log

For further assistance, please refer to the official documentation or contact support.
"""

    if use_tk:
        help_window = tk.Toplevel()
        help_window.title("SMATool Help")

        text_widget = tk.Text(help_window, wrap="word", width=100, height=30)
        text_widget.pack(expand=True, fill="both")

        text_widget.insert("1.0", help_text)
        text_widget.config(state="disabled")

        scrollbar = tk.Scrollbar(help_window, command=text_widget.yview)
        text_widget.config(yscrollcommand=scrollbar.set)
        scrollbar.pack(side="right", fill="y")
    else:
        print(help_text)


def smatoolguicall():
    root = tk.Tk()
    root.title("SMATool - Strength of Materials Analysis Toolkit - Input Generator")
    
    root.geometry("900x700")

    tabControl = ttk.Notebook(root)
    tab1 = ttk.Frame(tabControl)
    tabControl.add(tab1, text='Generate SMATool Input')
    tabControl.pack(expand=1, fill="both")

    home_dir = os.path.expanduser("~")
    smatool_dir = os.path.join(home_dir, "smatool_dir" if sys.platform == "win32" else ".smatool_dir")
    if not os.path.exists(smatool_dir):
        os.makedirs(smatool_dir)

    file_handler = FileHandler(smatool_dir)
    UFCFileGenerator(tab1, file_handler)

    root.mainloop()


if __name__ == "__main__":
    smatoolguicall()

