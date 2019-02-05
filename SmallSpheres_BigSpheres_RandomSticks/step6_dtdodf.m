
model         = 'dtd';
odfbin = 1; %sticks

% Prepare options
opt = mdm_opt();
opt = dtd_opt(opt);
opt.dtd.ind_start = 2;

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
NBS = 96;
odf_cell = cell(NBS,1);
for nBS = 1:NBS
    o     = fullfile(data_path,'NII_RES',model,'bootstrap',num2str(nBS));
    msf_mkdir(o);

    paths.nii_path = fullfile(o, 'maps');
    paths.mfs_fn   = fullfile(o, 'mfs.mat');
    paths.odf_fn   = fullfile(o, 'odf.mat');

    odf_cell{nBS} = dtd_4d_fit2odf(paths.mfs_fn, paths.odf_fn, opt);
    %mdm_fit2param(@dtd_4d_fit2odf, paths.mfs_fn, paths.odf_fn, opt);
end
%%
odf = odf_cell{1};

for nBS = 2:NBS
    odf_temp =  odf_cell{nBS};
    odf.w = cat(5,odf.w,odf_temp.w);
    for nbin = 1:numel(odf.w_bin)
        odf.w_bin{nbin} = cat(5,odf.w_bin{nbin},odf_temp.w_bin{nbin});
    end
end

odf_bsmean = odf;
odf_bsmean.w = mean(odf.w,5);
for nbin = 1:numel(odf.w_bin)
    odf_bsmean.w_bin{nbin} = mean(odf.w_bin{nbin},5);
end

%%
odf = odf_cell{1};

odf_s.n = odf.n;
odf_s.x = odf.x;
odf_s.y = odf.y;
odf_s.z = odf.z;
odf_s.c = abs([odf_s.x odf_s.y odf_s.z]);
odf_s.tri = odf.tri;

odf_s.w = squeeze(odf_bsmean.w_bin{odfbin}(9,10,1,:));
odf_s.verts = repmat(odf_s.w,[1 3]).*[odf_s.x odf_s.y odf_s.z];
odf_s.norms = vertexNormal(triangulation(odf_s.tri,odf_s.verts),(1:odf_s.n)');

figure(1), clf
p = patch('Faces',odf_s.tri,'Vertices',odf_s.verts);
axis tight, axis square, axis equal
view(0,0)
xlabel('x')
set(p,'FaceColor','interp','FaceVertexCData',odf_s.c,...
'EdgeColor','none','LineWidth',.1)
%axis off


sz = size(odf_bsmean.w);

odf_1 = odf_s;

odf_array.verts = [];
odf_array.norms = [];
odf_array.tri = [];
odf_array.c = [];

odf_array.wmax = max(odf_bsmean.w_bin{odfbin}(:));

wnorm = 2*odf_array.wmax;
dx = 1;
dy = dx;
dz = dx;

count = 0;
for nk = 1:sz(3)
    for nj = 1:sz(2)
        for ni = 1:sz(1)
%     for nj = 7:14
%         for ni = 7:14

            odf_2.n = odf.n;
            odf_2.x = odf.x;
            odf_2.y = odf.y;
            odf_2.z = odf.z;
            odf_2.c = abs([odf_s.x odf_s.y odf_s.z]);
            odf_2.tri = odf.tri;

            odf_2.w = squeeze(odf_bsmean.w_bin{odfbin}(ni,nj,nk,:))/wnorm;
            
            if sum(odf_2.w)>0
                odf_2.verts = repmat(odf_2.w,[1 3]).*[odf_2.x odf_2.y odf_2.z];
                odf_2.verts = odf_2.verts + [ni*ones(odf.n,1)*dx nj*ones(odf.n,1)*dy nk*ones(odf.n,1)*dz];
                odf_2.norms = vertexNormal(triangulation(odf_2.tri,odf_2.verts),(1:odf_2.n)');

                odf_array.verts = cat(1,odf_array.verts,odf_2.verts);
                odf_array.tri = cat(1,odf_array.tri,odf_2.tri+count*odf.n);
                odf_array.c = cat(1,odf_array.c,odf_2.c);
                odf_array.norms = cat(1,odf_array.norms,odf_2.norms);

                count = count+1;
            end
        end
    end
end

figure(2), clf
p = patch('Faces',odf_array.tri,'Vertices',odf_array.verts);
axis tight, axis square, axis equal
view(-20,30)
xlabel('x')
set(p,'FaceColor','interp','FaceVertexCData',odf_array.c,...
'EdgeColor','none','LineWidth',.1)
%axis off


odf_path = fullfile(data_path,'NII_RES',model,'odf');
msf_mkdir(odf_path);

fid = fopen(fullfile(odf_path,'ODF_verts.txt'), 'w');
N = size(odf_array.verts, 1);
fprintf(fid, '%8.0i,\n', N);
format = '%8.5f, %8.5f, %8.5f,\n';
for n = 1:N
    fprintf(fid,format,odf_array.verts(n,:));
end
fclose(fid);

fid = fopen(fullfile(odf_path,'ODF_norms.txt'), 'w');
N = size(odf_array.norms, 1);
fprintf(fid, '%8.0i,\n', N);
format = '%8.5f, %8.5f, %8.5f,\n';
for n = 1:N
    fprintf(fid,format,odf_array.norms(n,:));
end
fclose(fid);

fid = fopen(fullfile(odf_path,'ODF_c.txt'), 'w');
N = size(odf_array.c, 1);
fprintf(fid, '%8.0i,\n', N);
format = '%8.5f, %8.5f, %8.5f,\n';
for n = 1:N
    fprintf(fid,format,odf_array.c(n,:));
end
fclose(fid);

fid = fopen(fullfile(odf_path,'ODF_tri.txt'), 'w');
N = size(odf_array.tri, 1);
fprintf(fid, '%8.0i,\n', N);
format = '%8.0f, %8.0f, %8.0f,\n';
for n = 1:N
    fprintf(fid,format,odf_array.tri(n,:)-1);
end
fclose(fid);


toc;



