% Driver for analyzing DT_axderare2d data
% Execute from the opt/data/<user>/nmr/<dataset>/<expno> folder
% containing the ser and acqus files

% Select models
% 1) dti_euler
%    Conventional per-voxel diffusion tensors
%
% 2) dtd_gamma
%    Mean-square anisotropy and variance of isotropic diffusivities
%    Lasic et al, Front. Phys. 2, 11 (2014) 
%    http://dx.doi.org/10.3389/fphy.2014.00011
%   
% 3) dtd_covariance
%    Conventional per-voxel diffusion tensors
%    Mean-square anisotropy and variance of isotropic diffusivities
%    Westin et al, Neuroimage 135, 345-362 (2016) 
%    http://dx.doi.org/10.1016/j.neuroimage.2016.02.039
%   
% 4) dtd
%    Size-shape-orientation diffusion tensor distributions
%    Topgaard. NMR Biomed. e4066, (2019)
%    https://doi.org/10.1002/nbm.4066

clear

models         = {'dti_euler', 'dtd_gamma', 'dtd_covariance', 'dtd'};
c_model        = 1:4;

% Prepare paths
data_path = pwd;
in_path     = fullfile(data_path, 'indata'); 
out_path     = fullfile(data_path, 'fitdata');

msf_mkdir(out_path);

% Prepare options
opt = mdm_opt();
opt.verbose       = 1;
opt.do_overwrite  = 1;
opt.do_mask = 1; 
opt.mask.do_overwrite = 1;
opt.mask.threshold = .1;
opt.mask.b0_ind = 1;    
opt.filter_sigma = 0;
opt.do_data2fit = 1; 
opt.do_fit2param = 1;

% Connect to data
s.nii_fn = fullfile(in_path, 'data.nii.gz');
s.mask_fn = fullfile(in_path, 'data_mask.nii.gz');
s.xps = mdm_xps_load(fullfile(in_path, 'data_xps.mat'));

% Run analysis
for n_model = 1:numel(c_model)
    tic;
    
    % OUTPUT: define paths for data, fit parameters, and maps
    paths.nii_path = fullfile(out_path, 'maps');
    paths.mfs_fn   = fullfile(out_path, models{c_model(n_model)}, 'mfs.mat');
    paths.dps_fn   = fullfile(out_path, models{c_model(n_model)}, 'dps.mat');
    
    switch (models{c_model(n_model)})
        case 'dti_euler'
            nii_fn = dti_euler_pipe(s, paths, opt);
        case 'dtd_gamma'
            nii_fn = dtd_gamma_pipe(s, paths, opt);
        case 'dtd_covariance'
            nii_fn = dtd_covariance_pipe(s, paths, opt);
        case 'dtd'
            nii_fn = dtd_pipe(s, paths, opt);
    end
    
    toc;
end



