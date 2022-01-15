# Soundiation-Acoustophoresis
Soundiation: a MATLAB GUI used to predict the acoustic radiation force and torque, thereby the acoustophoresis of any axisymmetric particle under a plane wavefield or a user-specified transducer array.

Major features:
- The particle can be arbitrarily axisymmetric geometry;
- The sound-hard (Neumann) and sound-soft (Dirichlet) boundary conditions are provided;
- The medium can be air or water;
- The arrangment of transducer array, as well as phase & amplitude transducer parameters support user-specified (optional).

Theoretical background:
- The partial wave expansion method [1];
- The translation addiation theorem [2, 3];
- The conformal mapping technique [4, 5];

## Initialization

The program support parallel computation.

The source codes are found in the ```src``` folder. To start the Soundiation GUI, add the folder path:

``` matlab
addpath('<folderpath>/Soundiation-Acoustophoresis-main/src');
```

Then, open the the Soundiation GUI by typing in Command Window (MATLAB):

``` matlab
main_interface;
```

## Requirements

Soundiation has been tested with MATLAB2010a and above and should run on most personal laptops and desktop machines.

## Documentation

- ```./docs``` folder contains:
  -  an user manual;
  -  a code description.


## Functionality

Major functionality includes:
- Prediction of the acoustic radiation force and torque for non-spherical particles;
- Prediction of the translational and rotational motions of non-spherical particles above an user-specified transducer array.


## Contact
Tianquan Tang

tianquan@connect.hku.hk

## References

[1] E. G. Williams, Fourier acoustics: sound radiation and nearfeld acoustical holography, Academic Press, 1999, Chapter 6.

[2] P. A. Martin, Multiple scattering: interaction of time-harmonic waves with N obstacles, Cambridge University Press, 2006.

[3] T. Tang, L. Huang, Acoustic radiation force for multiple particles over a wide size-scale by multiple ultrasound sources, Journal of Sound and Vibration (2021) 116256.

[4] T. Tang, L. Huang, A fast semi-analytical procedure to calculate acoustic radiation force and torque for axisymmetric irregular bodies, xxx (2022).

[5] T. Tang, L. Huang, An efficient procedure to predict the acoustophoresis of axisymmetric irregular particles above ultrasound transducer array, xxx (2022).
