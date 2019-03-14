% Reconstruct nii from Bruker format
% Execute from the opt/data/<user>/nmr/<dataset>/<expno>/pdata_mddmri folder

% Image recon parameters
rps.smooth          = 300e-6; %Gaussian smoothing in [m]
rps.npix.read       = 16; % Number of pixels in read dimension
rps.npix.phase      = rps.npix.read; % Number of pixels in phase dimension
rps.shift_read      = -.25e-3; % Image shift in read dimension [m]

%-------------------------------------------------
wd = pwd;
data_path = fileparts(wd);
out_path  = fullfile(wd, 'indata'); 
msf_mkdir(out_path);

mdm_bruker_dt_rare2d_recon(data_path, out_path, rps);
