"""
  THICK-2D -- Thickness Hierarchy Inference & Calculation Kit for 2D materials

  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation
  version 3 of the License.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE.  See the GNU General Public License for more details.
  Email: cekuma1@gmail.com

""" 

import os

class WildCard:
    def __init__(self, settings_dict):
        self.settings_dict = settings_dict
        self.backup_dict = settings_dict.copy()

    def modify(self, key, value):
        if key in self.settings_dict:
            self.settings_dict[key] = value
        else:
            raise KeyError(f"{key} not found!")

    def reset(self):
        self.settings_dict.update(self.backup_dict)
        self._write_to_incar(self.settings_dict)

    def _write_to_incar(self, settings):
        with open("INCAR", "w") as f:
            for key, value in settings.items():
                f.write(f"{key.upper()} = {value}\n")



class IncarModifier:
    def __init__(self):
        self.original_values = {}

    def modify_incar(self, params=None, reset=False):
        #print("modify_incar called")  # Debugging
        
        if not params and not reset:
            return
        
        with open('INCAR', 'r') as f:
            lines = f.readlines()
        
        #print(f"Original INCAR: {lines}")  # Debugging
        existing_keys = set()
        
        for i, line in enumerate(lines):
            key = line.split('=')[0].strip()
            existing_keys.add(key)
            
            if reset:
                if key in self.original_values:
                    lines[i] = f' {key} = {self.original_values[key]}\n'
            else:
                if key in params:
                    if key not in self.original_values:
                        self.original_values[key] = line.split('=')[1].strip()
                        #print(f"Modifying key {key} to value {params[key]}")  # Debugging
                    lines[i] = f' {key} = {params[key]}\n'
        
        #print(f"Modified INCAR: {lines}")  # Debugging
        
        if not reset:
            for key, value in params.items():
                if key not in existing_keys:
                    lines.append(f' {key} = {value}\n')
        
        with open('INCAR', 'w') as f:
            f.writelines(lines)


class ChangeDir:
    def __init__(self, path):
        self.path = path
        self.original_path = os.getcwd()

    def __enter__(self):
        os.chdir(self.path)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.original_path)
