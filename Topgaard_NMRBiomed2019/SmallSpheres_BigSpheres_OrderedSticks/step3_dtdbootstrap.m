
model         = 'dtd';

% Prepare options
opt = mdm_opt();
opt = dtd_opt(opt);
opt.dtd.ind_start = 1;

opt.do_mask      = 0;
opt.do_data2fit  = 1;
opt.do_bootstrap = 1;
opt.do_datafit2chisq = 1;
opt.do_fit2param = 1;
opt.do_param2nii = 0;

opt.do_xps2pdf    = 0;
opt.do_nii2pdf    = 0;
opt.do_m2pdf      = 0;

opt.verbose       = 1;
opt.do_overwrite  = 1;

opt.dtd.do_plot = 0;
opt.mio.no_parfor = 1;

% Prepare paths
data_path = pwd;
i     = fullfile(data_path, 'NII_XPS'); 

% Connect to data
s.nii_fn = fullfile(i, 'data.nii.gz');
s.mask_fn = fullfile(i, 'data_mask.nii.gz');
s.xps = mdm_xps_load(fullfile(i, 'data_xps.mat'));


% Run analysis
tic;

% OUTPUT: define paths for data, fit parameters, and maps
parfor nBS = 1:96
    o     = fullfile(data_path,'NII_RES',model,'bootstrap',num2str(nBS));
    msf_mkdir(o);

    paths.nii_path = fullfile(o, 'maps');
    paths.mfs_fn   = fullfile(o, 'mfs.mat');
    paths.chisq_fn   = fullfile(o, 'chisq.mat');
    paths.ind_fn = fullfile(o, 'ind.mat');
    paths.dps_fn   = fullfile(o, 'dps.mat');
    paths = mdm_paths(paths);

    nii_fn = dtd_pipe(s, paths, opt);          
end

toc;



