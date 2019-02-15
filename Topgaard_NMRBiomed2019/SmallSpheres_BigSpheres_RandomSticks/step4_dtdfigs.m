% Generate figures corresponding to Fig 3, 5, and 6 in 
% Topgaard, Diffusion tensor distribution imaging, NMR Biomed. (2019)
% https://doi.org/10.1002/nbm.4066

%---------------------------------------------------------------------
clear all

% Parameter limits
plim.mdiso = 1.2e-9*[0 1];
plim.msqddelta = 1*[0 1];
plim.vdiso = 2.5e-19*[0 1];
plim.vsqddelta = .15*[0 1];
plim.cvdisosqddelta = 1.5e-10*[-1 1];

% slice
nk = 1;
% pixels
pix_i = [9 9 9];
pix_j = [10 7 4]; 
pix_k = nk*[1 1 1];

% pixel color
pix_col = [0 0.8 0
         1 .3 .3
         0 0 1];

% Font sizes etc
figscale = 2;
figwidth = figscale*17.78;
fs = figscale*6;
lw = figscale*1;
aspect = 1.618;
ms_max = figscale*20;
ms = figscale*3;

% Setup paths
model = 'dtd';
opt = mdm_opt();
opt = dtd_opt(opt);

data_path = pwd;
nii_xps_path = fullfile(data_path, 'NII_XPS');
nii_fn = fullfile(nii_xps_path, 'data.nii.gz');
xps_fn = fullfile(nii_xps_path, 'data_xps.mat');

res_path = fullfile(data_path, 'NII_RES');
paths = mdm_paths(res_path);
paths.mfs_fn = fullfile(fileparts(paths.mfs_fn), model, 'mfs.mat');
paths.dps_fn = fullfile(fileparts(paths.mfs_fn), 'dps.mat');
paths.maps = fullfile(paths.nii_path, 'maps');
bs_path = fullfile(res_path, model, 'bootstrap');

% Read data
[I,h] = mdm_nii_read(nii_fn);
xps = mdm_xps_load(xps_fn);
mfs = mdm_mfs_load(paths.mfs_fn);
dps = mdm_dps_load(paths.dps_fn);

%% Figure 3
% S0 with labeled pixels
figure(3), clf

im2d_S0 = dps.s0(:,:,nk);
s0max = max(reshape(dps.s0,numel(dps.s0),1));

height = .27;
width = height/aspect;
axh_S0_pixel = axes('position',[1-width 1-height-.12 width height]);
imagesc(im2d_S0'), hold on
colormap('gray')

set(axh_S0_pixel,'YDir','normal')
axis(axh_S0_pixel,'equal','tight','off')

for axh_n = 1:numel(axh_S0_pixel)
    for pix_n = 1:numel(pix_i)
        h = plot(axh_S0_pixel(axh_n),pix_i(pix_n),pix_j(pix_n),'o');
        set(h,'Color',pix_col(pix_n,:),'LineWidth',1*lw,'MarkerSize',ms)
    end
end


% Acquisition protocol
height = .1;
width = .7;
left = .05;
bottom_v = .52 + 1.14*[0*height 1*height 2*height 3*height];

xps_array = [round(xps.b,-7), round(xps.b_delta,1), xps.theta/pi*180, xps.phi/pi*180];
[xps_array_sort, ind_sort] = sortrows(xps_array,[1 2 3 4]);

axh_b = axes('position',[left bottom_v(4) width height]);
ph_b = plot(axh_b,1:xps.n,xps.b(ind_sort)/1e9);
axh_bdelta = axes('position',[left bottom_v(3) width height]);
ph_bdelta = plot(axh_bdelta,1:xps.n,xps.b_delta(ind_sort));
axh_theta = axes('position',[left bottom_v(2) width height]);
ph_theta = plot(axh_theta,1:xps.n,xps.theta(ind_sort)/pi*180);
axh_phi = axes('position',[left bottom_v(1) width height]);
ph_phi = plot(axh_phi,1:xps.n,xps.phi(ind_sort)/pi*180);

axh_prot_v = [axh_b; axh_bdelta; axh_theta; axh_phi];
ph_prot_v = [ph_b; ph_bdelta; ph_theta; ph_phi];

% Signals and DTDs for selected pixels
w_threshold = .01;
s0_threshold = .1;

height = .12;
left_v = left + [0 0 0];
bottom_v = .06 + 1.15*[2*height height 0];

dmin = opt.dtd.dmin;
dmax = opt.dtd.dmax;
ratiomin = dmin/dmax;
ratiomax = dmax/dmin;

xmin = log10(dmin);
xmax = log10(dmax);
ymin = log10(1/ratiomax);
ymax = log10(ratiomax);

axh_pix_v = [];
axh_dist_v = [];
for pix_n = 1:numel(pix_i)
    ni = pix_i(pix_n);
    nj = pix_j(pix_n);
    left = left_v(pix_n);
    bottom = bottom_v(pix_n);

    s0 = im2d_S0(ni,nj,nk);

    signal = squeeze(I(ni,nj,nk,:))/s0;
    signal_fit = dtd_1d_fit2data(squeeze(mfs.m(ni,nj,nk,:))', xps)/s0;
    rmsd = sqrt(sum((signal_fit-signal).^2)/numel(signal));
    
    axh_pix = axes('position',[left bottom width height]);
    h = plot(1:numel(signal),signal(ind_sort)','o');
    set(h,'Color',pix_col(pix_n,:),'LineWidth',1*lw,'MarkerSize',ms)
    hold on
    plot(1:numel(signal),signal_fit(ind_sort)','k.','MarkerSize',1.2*ms)
    axh_pix_v = [axh_pix_v; axh_pix];

    m = squeeze(mfs.m(ni,nj,nk,:))';
    s0 = dps.s0(ni,nj,nk);
    dtd = dtd_m2dtd(m);
    [n,par,perp,theta,phi,w] = dtd_dist2par(dtd);
    if n>0
        xcos = cos(phi).*sin(theta);
        ycos = sin(phi).*sin(theta);
        zcos = cos(theta);

        iso = tm_eigvals2iso([par perp perp]);
        fa = tm_eigvals2fa([par perp perp]);

        c.x = log10(iso);
        c.y = log10(par./perp);
        c.ms = ms_max*sqrt(w/s0max);
        c.bright = fa;
        c.r = abs(xcos);
        c.g = abs(ycos);
        c.b = abs(zcos);

        axh_dist = axes('position',[1-height/aspect-.04 bottom height/aspect height]);
        for nc = 1:n
            if w(nc) > w_threshold*s0
                h1 = plot(c.x(nc),c.y(nc),'o');
                hold on
                col = c.bright(nc)*[c.r(nc) c.g(nc) c.b(nc)];
                %set(h1,'MarkerSize',c.ms(nc),'Color',col,'MarkerFaceColor',col)
                set(h1,'MarkerSize',c.ms(nc),'Color',col,'MarkerFaceColor','none','LineWidth',.75*lw)
            end
        end
        axh_dist_v = [axh_dist_v; axh_dist];
    end
end

set(axh_b,'YLim',max(xps.b/1e9)*[-.1 1.1])
set(axh_bdelta,'YLim',[-.65 1.3],'YTick',[-.5 0 .5 1])
set(axh_theta,'YLim',180*[-.1 1.1],'YTick',[0:90:180])
set(axh_phi,'YLim',180*[-1.2 1.2],'YTick',[-180:90:180])
set(ph_prot_v,'Marker','.','LineStyle','none','Color','k','MarkerSize',1.2*ms)

set([axh_prot_v; axh_pix_v],'Box','off','TickDir','out','XLim',xps.n*[-.01 1.01],...
    'XTick',[0:100:500],'XTickLabel',[],'TickLength',.005*[1 1],'LineWidth',lw,'FontSize',fs)
set(axh_pix,'XTickLabel',[0:100:500])
set(axh_pix_v,'YLim',[-.05 1.05],'YTick',0:.5:1)

set(axh_dist_v,'XLim',[xmin xmax]+.01*(xmax-xmin)*[-1 1], 'YLim',[ymin ymax]+.01*(ymax-ymin)*[-1 1],'XTick',-11:-8,'YTick',-2:2,...
    'TickDir','out','TickLength',.025*[1 1],'LineWidth',lw,'FontSize',fs,'YAxisLocation','right','XTickLabel',[])
set(axh_dist,'XTickLabel',-11:-8)

set(gcf, 'PaperUnits','centimeters','PaperPosition', figwidth*[0 0 1 1/aspect],'PaperSize', figwidth*[1 1/aspect]);
fig_fn = fullfile(paths.maps,[model '_protocol_pixelssignals']);
eval(['print ' fig_fn ' -dpdf -loose'])

%% Figure 5
param_dist = {'dpar','dperp','theta','phi','w'};
param_map = {'s0','mdiso','msqddelta','vdiso','vsqddelta','cvdisosqddelta'};
param = {param_dist{:}, param_map{:}, 'chisq'};
for nparam = 1:numel(param)
    eval([param{nparam} ' = [];']);
end
for nparam = 1:numel(param_map)
    eval([param_map{nparam} ' = [];']);
end
for nbin = 1:numel(opt.dtd.bin_disomax)
    bin{nbin}.no = nbin;
    for nparam = 1:numel(param_map)
        eval(['bin{nbin}.' param_map{nparam} ' = [];']);
    end
end

Nbins = numel(bin);
Nparams = numel(param_map);

bsno = msf_getdirno(bs_path);
nn = 0;
for nbs = bsno
    mfs_fn   = fullfile(bs_path, num2str(nbs), 'mfs.mat');
    dps_fn   = fullfile(bs_path, num2str(nbs), 'dps.mat');
    chisq_fn   = fullfile(bs_path, num2str(nbs), 'chisq.mat');
    if exist(mfs_fn,'file')==2
        temp = load(mfs_fn); mfs = temp.mfs; 
        m = mfs.m;
        sz = size(m);
        n = m(:,:,:,1);
        nn_temp = (sz(4)-1)/6;

        ind = false(sz(4),1);
        ind(2:5:end) = 1;

        for nparam = 1:numel(param_dist)
            eval([param_dist{nparam} ' = cat(4,' param_dist{nparam} ',m(:,:,:,circshift(ind,' num2str(nparam-1) ',1)));']);
        end
        nn = nn + nn_temp;
    end
    if exist(dps_fn,'file')==2
        temp = load(dps_fn); dps = temp.dps; 
        for nparam = 1:numel(param_map)
            eval([param_map{nparam} ' = cat(4,' param_map{nparam} ',dps.' param_map{nparam} ');']);
        end
        for nbin = 1:numel(dps.bin)
            for nparam = 1:numel(param_map)
                eval(['bin{nbin}.' param_map{nparam} ' = cat(4,bin{nbin}.' param_map{nparam} ',dps.bin{nbin}.' param_map{nparam} ');']);
            end
        end
    end
    if exist(chisq_fn,'file')==2
        temp = load(chisq_fn); 
        chisq = cat(4,chisq,temp.chisq);
    end
end

for nparam = 1:numel(param_map)
    eval(['dps_bsmean.' param_map{nparam} ' = mean(' param_map{nparam} ', 4);']);
    for nbin = 1:numel(dps.bin)
        eval(['dps_bsmean.bin{nbin}.' param_map{nparam} ' = mean(bin{nbin}.' param_map{nparam} ', 4);']);
    end
end

figure(5), clf
axh_v = dtd_dpsbins2maps(gcf,dps_bsmean,plim);
papersize = 3*[Nparams Nbins+1];
set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
fig_fn = fullfile(paths.maps,[model '_pmaps']);
eval(['print ' fig_fn ' -dpdf -loose'])

%% Figure 6
nbin_tot = numel(bin)+1;
for nparam = 1:numel(param_map)
    eval(['bin{nbin_tot}.' param_map{nparam} ' = ' param_map{nparam} ';']);
end

bin = bin(fliplr(circshift(1:nbin_tot,[1,1])));

s0max = max(dps_bsmean.s0(:));
s0thresh = 0.01;

left = .05; bottom = .07; dleft = (1-left)/Nparams; dheight = (1-bottom)/(Nbins+1);
width = .8*dleft; height = .9*dheight;

Nhistbins = 100;

plim.s0 = s0max*[0 1];
pscale.s0 = s0max;
pscale.mdiso = 1e-9;
pscale.msqddelta = 1;
pscale.vdiso = 1e-18;
pscale.vsqddelta = 1;
pscale.cvdisosqddelta = 1e-9;

figure(6), clf

for nparam = 1:Nparams
    for nbin = 1:(Nbins+1)
        dps_temp = bin{nbin};
        eval(['mapdat = dps_temp.' param_map{nparam} '/pscale.' param_map{nparam} ';']);
        eval(['xlim = plim.' param_map{nparam} '/pscale.' param_map{nparam} ';']);
        axh = axes('position',[left+(nparam-1)*dleft (nbin-1)*dheight+bottom width height]);
        edges = linspace(min(xlim),max(xlim),Nhistbins);
        for npix = 1:3
            ni = pix_i(npix); nj = pix_j(npix); nk = pix_k(npix);
            col = pix_col(npix,:);
            histdat = squeeze(mapdat(ni,nj,nk,:));
            histdats0 = squeeze(dps_temp.s0(ni,nj,nk,:));
            if ~strcmp(param_map{nparam},'s0')
                histdat = histdat(histdats0>s0thresh*s0max);
            end
            h = histogram(histdat,edges,'DisplayStyle','stairs');
            set(h,'EdgeColor',col,'LineWidth',.25*lw*(4-npix))
            set(h,'LineWidth',lw)
            hold on
        end
        set(axh,'XLim',xlim,'YLim',100*[-.1 1.1],'Box','off','TickDir','out','TickLength',.03*[1 1],'LineWidth',lw,'FontSize',fs)
        if nbin > 1
            set(axh,'XTickLabel',[])
        end
        if nparam > 1
            set(axh,'YTickLabel',[])
        end
    end
end

papersize = figwidth*[1 1/aspect];
set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
fig_fn = fullfile(paths.maps,[model '_pixelshistograms']);
eval(['print ' fig_fn ' -dpdf -loose'])

 