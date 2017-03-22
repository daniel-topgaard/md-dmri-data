# Multidimensional diffusion MRI data
This repository contains example data that can be processed with the [MD-dMRI framework](https://github.com/daniel-topgaard/md-dmri/). The data was acquired on a Bruker microimaging system using tubes with liquid crystals as described in [Topgaard 2016](http://dx.doi.org/10.1039/c5cp07251d) and [Martins 2016](http://dx.doi.org/10.1103/PhysRevLett.116.087601). The following example folders are included:
* [OrderedSticks](OrderedSticks). See fig 4 (top right) in [Topgaard 2016](http://dx.doi.org/10.1039/c5cp07251d).
* [RandomSticks](RandomSticks). See fig 4 (bottom right) in [Topgaard 2016](http://dx.doi.org/10.1039/c5cp07251d).
* [OrderedPlanes](OrderedPlanes). See fig 4 (top left) in [Topgaard 2016](http://dx.doi.org/10.1039/c5cp07251d).
* [SmallSpheres_BigSpheres_RandomSticks](SmallSpheres_BigSpheres_RandomSticks). Sample geometry is shown in fig 3 in [Martins 2016](http://dx.doi.org/10.1103/PhysRevLett.116.087601). 
* [SmallSpheres_BigSpheres_RandomSticks](SmallSpheres_BigSpheres_RandomSticks). As above, but with magnetically aligned liquid crystalline domains. 

## How to start
Run setup_paths.m in the root folder to put files in the Matlab path. Run step1_recon_data.m and step2_run_analysis.m in the example folders above to process the data. View the raw or processed data with the GUI included in the [MD-dMRI framework](https://github.com/daniel-topgaard/md-dmri/). The GUI is started by typing mgui in the Matlab command window.
