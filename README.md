# Multidimensional diffusion MRI data
This repository contains some simple data sets suitable for getting acquainted with the [MD-dMRI framework](https://github.com/daniel-topgaard/md-dmri/) for processing multidimensional diffusion MRI data.<sup>1</sup> The data was acquired on a Bruker microimaging system with a custom diffusion-encoded RARE sequence<sup>2</sup> using gradient waveforms and q-trajectories<sup>3</sup> giving axisymmetric b-tensors.<sup>4</sup> The data was originally reported in
* [Topgaard_NMRBiomed2019](https://doi.org/10.1002/nbm.4066)
* [Jiang_MScProject2018](http://www.physchem.lu.se/people/phdstudents/jiang/) 
* [Topgaard_PhysChemChemPhys2016](http://dx.doi.org/10.1039/c5cp07251d).

## How to start
* Download the [MD-dMRI framework](https://github.com/daniel-topgaard/md-dmri/) and add to the Matlab search path.
* Download this repository and locate a pdata_mddmri folder within for instance the Topgaard_NMRBiomed2019/SmallSpheres_BigSpheres_OrderedSticks data set.
* Run step1_recon.m to convert the raw Bruker data to nifti images and experimental parameters to matlab format.
* Run step2_fit.m to process the data with gamma,<sup>2</sup> covariance,<sup>5</sup> and 4D diffusion tensor distribution.<sup>6</sup>
* View the images in the indata and fitdata/maps folders with the GUI started by typing `mgui` in the Matlab command window.
* Run step3_dtdbootstrap.m to perform uncertainty estimation for the 4D diffusion tensor distributions.<sup>6</sup>
* Run step4_dtdfigs.m to generate figures corresponding to Figs 3-7 in [Topgaard_NMRBiomed2019](https://doi.org/10.1002/nbm.4066).

# References
1. D. Topgaard. Multidimensional diffusion MRI. [J. Magn. Reson. 275, 98-113 (2017)](http://dx.doi.org/10.1016/j.jmr.2016.12.007).
2. S. Lasič, F. Szczepankiewicz, S. Eriksson, M. Nilsson, D. Topgaard. Microanisotropy imaging: quantification of microscopic diffusion anisotropy and orientational order parameter by diffusion MRI with magic-angle spinning of the q-vector. [Front. Physics 2, 11 (2014)](http://dx.doi.org/10.3389/fphy.2014.00011).
3. S. Eriksson, S. Lasič, D. Topgaard. Isotropic diffusion weighting by magic-angle spinning of the q-vector in PGSE NMR. [J. Magn. Reson. 226, 13-18 (2013)](http://dx.doi.org/10.1016/j.jmr.2012.10.015).
4. S. Eriksson, S. Lasič, M. Nilsson, C.-F. Westin, D. Topgaard. NMR diffusion encoding with axial symmetry and variable anisotropy: Distinguishing between prolate and oblate microscopic diffusion tensors with unknown orientation distribution. [J. Chem. Phys. 142, 104201 (2015)](http://dx.doi.org/10.1063/1.4913502).
5. C.-F. Westin, H. Knutsson, O. Pasternak, F. Szczepankiewicz, E. Özarslan, D. van Westen, C. Mattisson, M. Bogren, L. O'Donnell, M. Kubicki, D. Topgaard, M. Nilsson. Q-space trajectory imaging for multidimensional diffusion MRI of the human brain. [Neuroimage 135, 345-362 (2016)](http://dx.doi.org/10.1016/j.neuroimage.2016.02.039).
6. D. Topgaard. Diffusion tensor distribution imaging. [NMR Biomed., e4066 (2019)](https://doi.org/10.1002/nbm.4066).