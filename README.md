# Soundiation-Radiation-force-torque-Acoustophoresis
A MATLAB GUI, Soundiation, used to predict the acoustic radiation force and torque, thereby the acoustophoresis of any axisymmetric particle under a plane wavefield or a user-specified transducer array.

Major features:
- The particle can be arbitrarily axisymmetric geometry;
- The sound-hard (Neumann) and sound-soft (Dirichlet) boundary conditions are provided;
- The medium can be air or water;
- The arrangment of transducer array, as well as phase & amplitude transducer parameters support user-specified (optional).

Theoretical background:
- The partial wave expansion method [];
- The translation addiation theorem [];
- The conformal mapping technique [];

## Initialization

The program support parallel computation.

To start the Soundiation GUI, add the folder path:

``` matlab
addpath('<folderpath>/Soundiation/src');
```

The source code is found in the ```src``` folder.

To open the the Soundiation GUI, type:

``` matlab
main_interface;
```

## Requirements

Soundiation has been tested with MATLAB2010a and above and should run on most personal laptops and desktop machines.

## Documentation

- To view the html documentation in MATLAB:
  - ```help (or press F1) -> Supplemental Software -> ElasticMatrix Toolbox```

- ```./examples``` folder contains example scripts demonstrating some of the capabilities of the code.
- ```./examples_mlx``` folder contains example scripts in the MATLAB live script style.
- ```./documentation``` folder contains a reference list for the mathematical background.

## Functionality

Major functionality includes:
- Plotting slowness profiles.


## Contact
Tianquan Tang

tianquan@connect.hku.hk

## References

[1] 
