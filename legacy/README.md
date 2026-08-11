# Proply (legacy Python version)

This is the original Python implementation of proply, superseded by the Rust
port in `proply-rs/` (see the [top-level README](../README.md)). It is kept
for reference: the BEM equations, NACA/ARA-D foil families, and the design
loop were ported 1:1 to Rust, and `build/pyref/` + the golden tests still
use this package as the reference implementation.

The Python workflow designs a blade as an STL mesh and generates an
OpenSCAD model; the Rust port writes a single STEP (AP242) file instead.

## Install

Proply depends on a slightly modified version of xfoil-python:

    git clone https://github.com/mxjeff/xfoil-python.git
    cd xfoil-python
    pip3 install .

Then install proply itself:

    make develop

## Creating a blade

Specify the prop with a JSON file (see `../props/test_prop.json` for an
example), then run:

    make test TARGET=test_prop

or directly:

    proply --naca --bem --n 40 --resolution 30 --param='../props/test_prop.json'

The `make test` target runs meshlabserver on the generated STL to clean up
duplicate vertices.

## Docker

    make build
    make run

See `Makefile` and `docker-compose.yml` for details.
