# Proply: Propeller Design in Python

Note: This software is a work-in-progress, its also under heavy development, so I may break things. Please accept my apologies if this happens.

Contact the Author (Tim Molteno) if you have any questions.

This code will generate propeller designs automatically that can be 3D printed. Below is an example of a three-bladed prop designed with this software.

![alt text][prop5x3]

## Install

Proply depends on a slightly modified version of xfoil-python. It should be installed as follows

    git clone https://github.com/mxjeff/xfoil-python.git
    cd xfoil-python
    pip3 install .

## Specifying a prop

This is done using a JSON file. Here is an example:

    {
        "name": "FLYWOO_ROBO_RB1202.5",
        
        "altitude": 0.0,
        "forward_airspeed": 0.0,
        
        "motor_Kv": 11500,
        "motor_volts": 3.7,
        "motor_no_load_current": 0.075,
        "motor_winding_resistance": 0.035,
        
        "thrust": 0.5,
        "blades": 3,

        "radius": 0.02,
        "tip_chord": 0.003,
        "center_hole": 1.5e-3,
        "scimitar_percent": -5.0,
        "trailing_edge": 0.25,
        "hub_radius": 0.004,
        "hub_depth": 0.003
    }

The prop will automatically be optimized for the motor to operate at max efficiency at the specified thrust. For this to work, you need to put in the motor parameters.

### Creating a blade

A single blade is generated as an STL file using the command:

    proply --naca --bem --n 40 --resolution 30 --param='test_prop.json'

The mesh can be cleaned using meshlab to remove duplicate vertices with:

    make scad TARGET=test_prop

##  Using OpenSCAD to generate a propeller

Proply generates an OpenSCAD file. This file imports the blade STL and generates a propeller model.

## TODO

* Optimize foils for angle of attack and shape.

# Creating Props using JavaProp, FreeCAD and OpenSCAD.

This is an alternative approach, allowing one to generate STL files from javaprop designs.

## Step 1. Design your prop in JavaProp.

## Step 2. Export your design in IGES format

## Step 3. Create and STL fvrom IGES

This is done using FreeCAD.

* Open FreeCAD and Select Part Mode
* Import your IGES file
* Choose Part -> Shape Builder
* Choose "Face From Edges" and select the boundary of one end of the prop blade. Click "Create"
* Do the same for the other end of the prop blade.
* Choose "Shell from Faces" and select (holding down CTRL) the two ends, and the surface of the prop blade. Click "Create"

You will now have a design with several parts. One of them is called "Shell". Select this, and choose Part -> Convert to solid.

now delete all the other bits of your design (right click and delect the original import, as well as the other parts.

Then export this as an .stl file.


# Use OpenSCAD

Import the .stl file inside openscad. The file prop5x3.scad shows how this is done.


[prop5x3]: https://github.com/elec-otago/proply/blob/master/images/prop5x3.png "Three Bladed Prop"
