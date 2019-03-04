% Generate odf figure corresponding to Fig 7 in 
% Topgaard, Diffusion tensor distribution imaging, NMR Biomed. (2019)
% https://doi.org/10.1002/nbm.4066
clear all

dtdbin = 4; %1:sticks, 4:planes

% Font sizes etc
figscale = 2;
figwidth = figscale*17.78;
fs = figscale*6;
lw = figscale*1;
aspect = 1.618;

% Prepare options
opt = mdm_opt();
opt = dtd_opt(opt);

% Prepare paths
model         = 'dtd';
data_path = pwd;
maps_path = fullfile(data_path, 'NII_RES', 'maps');
bs_path = fullfile(data_path,'NII_RES',model,'bootstrap');

% Collect data from bootstraps
bsno = msf_getdirno(bs_path);
nn = 0;
for nbs = bsno
    mfs_fn   = fullfile(bs_path, num2str(nbs), 'mfs.mat');
    odf_cell{nbs} = dtd_4d_fit2odf(mfs_fn, opt);
end

odf = odf_cell{1};

for nbs = bsno(2:end)
    odf_temp =  odf_cell{nbs};
    odf.w = cat(5,odf.w,odf_temp.w);
    for nbin = 1:numel(odf.w_bin)
        odf.w_bin{nbin} = cat(5,odf.w_bin{nbin},odf_temp.w_bin{nbin});
    end
end

% Mean over bootstraps
odf_bsmean = odf;
odf_bsmean.w = mean(odf.w,5);
for nbin = 1:numel(odf.w_bin)
    odf_bsmean.w_bin{nbin} = mean(odf.w_bin{nbin},5);
end

%%
sz = size(odf_bsmean.w);

odf_array.verts = [];
odf_array.norms = [];
odf_array.tri = [];
odf_array.c = [];

odf_array.wmax = max(odf_bsmean.w_bin{dtdbin}(:));

wnorm = 1.5*odf_array.wmax;
dx = 1;
dy = dx;
dz = dx;

count = 0;
for nk = 1:sz(3)
    for nj = 1:sz(2)
        for ni = 1:sz(1)

            odf_voxel.n = odf.n;
            odf_voxel.x = odf.x;
            odf_voxel.y = odf.y;
            odf_voxel.z = odf.z;
            odf_voxel.c = abs([odf.x odf.y odf.z]);
            odf_voxel.tri = odf.tri;

            odf_voxel.w = squeeze(odf_bsmean.w_bin{dtdbin}(ni,nj,nk,:))/wnorm;
            
            if sum(odf_voxel.w)>0
                odf_voxel.verts = repmat(odf_voxel.w,[1 3]).*[odf_voxel.x odf_voxel.y odf_voxel.z];
                odf_voxel.verts = odf_voxel.verts + [ni*ones(odf.n,1)*dx nj*ones(odf.n,1)*dy nk*ones(odf.n,1)*dz];
                odf_voxel.norms = vertexNormal(triangulation(odf_voxel.tri,odf_voxel.verts),(1:odf_voxel.n)');

                odf_array.verts = cat(1,odf_array.verts,odf_voxel.verts);
                odf_array.tri = cat(1,odf_array.tri,odf_voxel.tri+count*odf.n);
                odf_array.c = cat(1,odf_array.c,odf_voxel.c);
                odf_array.norms = cat(1,odf_array.norms,odf_voxel.norms);

                count = count+1;
            end
        end
    end
end

figure(2), clf
p = patch('Faces',odf_array.tri,'Vertices',odf_array.verts);
axis tight, axis square, axis equal
view(-20,30)
xlabel('x'), ylabel('y')
set(p,'FaceColor','interp','FaceVertexCData',odf_array.c,...
'EdgeColor','none','LineWidth',.1)
set(gca,'FontSize',fs,'LineWidth',lw)
%axis off

papersize = figwidth*[1 1/aspect];
set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
fig_fn = fullfile(maps_path,[model '_odf']);
eval(['print ' fig_fn ' -dpng -loose'])


return

%Export to povray

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



