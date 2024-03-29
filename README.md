# Soundiation-Acoustophoresis
Soundiation: a MATLAB GUI-based software used to predict the acoustic radiation force and torque, thereby the acoustophoresis of any axisymmetric particle under a plane wavefield or a user-customized transducer array.

![image](https://github.com/Tountain/Images/blob/main/GUI.bmp)

Major features:
- The particle can be designed as arbitrarily axisymmetric geometry (by the mapping coefficients "c_n");
- The sound-hard (Neumann) and sound-soft (Dirichlet) boundary conditions are provided;
- The medium can be air or water;
- The arrangment of transducer array, as well as phase & amplitude transducer parameters support user-specified (optional).

Theoretical background:
- The partial wave expansion method [1];
- The translation addition theorem [2, 3];
- The conformal mapping technique [4, 5];

## Initialization

The program support parallel computation (Refer to user manual in ```./doc/user_manual``` for details). 

The source codes are found in the ```./src``` folder. To start the Soundiation GUI, add the folder path:

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

- File ```arXiv_2202.04526``` is a description preprint for the software.

- ```./docs``` folder contains:
  -  A numerical model (by COMSOL Multiphysics 5.5) to validate the calculation results, if needed;
  -  A user manual.

- ```./src``` folder contains:
  -  All source codes (".m") and a GUI framework ("main_interface.fig") for the software.

- ```./data file (example)``` folder contains: 
  -  An example of the user designed particle geometry ("particle_data.stl");
  -  An example of the predicted dynamic data (namely "Myfilename.txt" by default).

   **Note**: By default, the above data files are automatically saved in the ```./src``` folder.


## Functionality

Major functionalities includes:
- Design a non-spherical particle and output a "particle_data.stl" file;
- Prediction of the acoustic radiation force and torque on non-spherical particles;
- Prediction of the dynamics (translational and rotational motions) of non-spherical particles above an user-specified transducer array (the dynamic data is saved in a "Myfilename.txt" by default).


## Contact
Tianquan Tang

- Email address: tianquan@connect.hku.hk; ttqtianquan@gmail.com.

- Researchgate: https://www.researchgate.net/profile/Tianquan-Tang.


## Major references and paper

- References

[1] E. G. Williams, Fourier acoustics: sound radiation and nearfeld acoustical holography, Academic Press, 1999, Chapter 6.

[2] P. A. Martin, Multiple scattering: interaction of time-harmonic waves with N obstacles, Cambridge University Press, 2006.

[3] T. Tang, L. Huang, Acoustic radiation force for multiple particles over a wide size-scale by multiple ultrasound sources, Journal of Sound and Vibration 509 (2021) 116256. DOI: https://doi.org/10.1016/j.jsv.2021.116256

[4] T. Tang, L. Huang, An efficient semi-analytical procedure to calculate acoustic radiation force and torque for axisymmetric irregular bodies, Journal of Sound and Vibration 532 (2022) 117012. DOI: https://doi.org/10.1016/j.jsv.2022.117012

[5] T. Tang, L. Huang, Theoretical framework to predict the acoustophoresis of axisymmetric irregular objects above an ultrasound transducer array, Physical Review E 105 (2022) 055110 (or arXiv:2201.13042). DOI: https://doi.org/10.1103/PhysRevE.105.055110

- Paper

[1] T. Tang, L. Huang, Soundiation: A software in evaluation of acoustophoresis driven by radiation force and torque on axisymmetric objects, The Journal of the Acoustical Society of America 152.5 (2022) 2934. DOI: https://doi.org/10.1121/10.0015199

