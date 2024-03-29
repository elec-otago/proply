"""
    PROPLY create printable props in python
    Author: Tim Molteno, tim@elec.ac.nz
    Copyright (c) 2019-2021.

    License. GPLv3.
"""
from setuptools import setup, find_packages

with open("README.md") as f:
    readme = f.read()

setup(
    name="proply",
    version="0.1.4b7",
    description="Convert TART observation data to Measurement Sets",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="http://github.com/elec-otago/proply",
    author="Tim Molteno",
    author_email="tim@elec.ac.nz",
    license="GPLv3",
    install_requires=["scipy", "matplotlib", "pyaml", "numpy", "numpy-stl",
                      "xfoil @ https://github.com/mxjeff/xfoil-python/tarball/master"],
    package_data={
        "proply": [
            "sql/foil_simulator.sql",
            "foils/foil_simulator.sql",
        ]
    },
    include_package_data=True,
    test_suite="nose.collector",
    tests_require=["nose"],
    packages=["proply", "proply.sql", "proply.foils", "proply.templates"],
    scripts=["bin/proply"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Intended Audience :: Science/Research",
    ],
)
