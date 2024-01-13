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

setup(
    name='qepotential',
    version='1.0',
    scripts=['download_qe_pseudo.py'],
    author='Chinedu Ekuma',
    author_email="cekuma1@gmail.com",
    description="Utility for SMATool - Automated toolkit for computing zero and finite-temperature strength of materials",
    license="GNU GPL version 3",
    install_requires=[
        'selenium',
        'webdriver-manager'
    ],
    entry_points={
        'console_scripts': [
            'qepotential=download_qe_pseudo:main'
        ],
    },
    include_package_data=True,
    zip_safe=False,
)

