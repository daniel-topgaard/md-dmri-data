% Generate figures corresponding to Fig 3, 4, 5, 6, and 7 in 
% Topgaard, Diffusion tensor distribution imaging, NMR Biomed. (2019)
% https://doi.org/10.1002/nbm.4066

%---------------------------------------------------------------------
clear all

model = 'dtd';

% Voxels
voxels.i = [9 9 9];
voxels.j = [10 7 4]; 
nk = 1;
voxels.k = nk*[1 1 1];

odfbin = 1; %1:sticks, 4:planes

% Parameter limits
plim.mdiso = 1.2e-9*[0 1];
plim.msddelta = 1*[0 1];
plim.vdiso = 2.5e-19*[0 1];
plim.vsddelta = .15*[0 1];
plim.cvdisosddelta = 1.5e-10*[-1 1];

% Options
opt = mdm_opt();
opt = dtd_opt(opt);
opt.dtd.plim = plim;
opt.dtd.odfbin = odfbin; %1:sticks, 4:planes

% Setup paths
wd = pwd;
paths.indata_path = fullfile(wd, 'indata');
paths.nii_fn = fullfile(paths.indata_path, 'data.nii.gz');
paths.xps_fn = fullfile(paths.indata_path, 'data_xps.mat');

paths.fitdata_path = fullfile(wd, 'fitdata');
paths.mfs_fn = fullfile(paths.fitdata_path, model, 'mfs.mat');
paths.dps_fn = fullfile(paths.fitdata_path, model, 'dps.mat');
paths.bs_path = fullfile(paths.fitdata_path, model, 'bootstrap');
paths.figs = fullfile(wd, 'dtdfigs');


%% Figure 3
mplot_dtd_protocol_voxelsignal(voxels, paths, opt)
%% Figure 4  
mplot_dtd_array_global(paths, opt)
%% Figures 5 and 6
mplot_dtd_pmaps_hist(voxels, paths, opt)
%% Figure 7
mplot_dtd_odf(paths, opt)
 