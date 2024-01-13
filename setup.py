"""
  SMATool - Automated toolkit for computing zero and finite-temperature strength of materials

  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation
  version 3 of the License.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  cekuma1@gmail.com
""" 

try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup


try:
    from importlib.metadata import version  # Python 3.8+
except ImportError:
    from importlib_metadata import version  # Python <3.8

#from setuptools import setup

setup()

