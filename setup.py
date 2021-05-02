'''
    PROPLY create printable props in python
    Author: Tim Molteno, tim@elec.ac.nz
    Copyright (c) 2019-2021.

    License. GPLv3.
'''
from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

setup(name='proply',
    version='0.1.4b5',
    description='Convert TART observation data to Measurement Sets',
    long_description=readme,
    long_description_content_type="text/markdown",
    url='http://github.com/elec-otago/proply',
    author='Tim Molteno',
    author_email='tim@elec.ac.nz',
    license='GPLv3',
    install_requires=['scipy', 'matplotlib', 'pyaml', 'numpy', 'numpy-stl'],
    package_data={'proply': ['foil_simulator.sql',]},

    test_suite='nose.collector',
    tests_require=['nose'],
    packages=['proply'],
    scripts=['bin/proply'],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        "Intended Audience :: Science/Research"])
