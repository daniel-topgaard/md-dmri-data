% Reconstruct nii from Bruker format
% Execute from the opt/data/<user>/nmr/<dataset>/<expno> folder
% containing the ser and acqus files

% Image recon parameters
rps.smooth          = 900e-6; %Gaussian smoothing in [m]
rps.npix.read       = 16; % Number of pixels in read dimension
rps.npix.phase      = rps.npix.read; % Number of pixels in phase dimension
rps.shift_read      = -.25e-3; % Image shift in read dimension [m]

%-------------------------------------------------
data_path = pwd;
out_path  = fullfile(data_path, 'NII_XPS'); 

msf_mkdir(out_path);

mdm_bruker_dt_rare2d_recon(data_path, out_path, rps);
