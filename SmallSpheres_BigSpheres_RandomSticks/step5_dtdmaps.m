clear all

model = 'dtd';
opt = mdm_opt();
opt = dtd_opt(opt);

% Prepare paths
data_path = pwd;
in_path     = fullfile(data_path, 'NII_XPS');
res_path     = fullfile(data_path, 'NII_RES');

param_dist = {'dpar','dperp','theta','phi','w'};
nn = 0;
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

Nparams = numel(param_map);
Nbins = numel(opt.dtd.bin_disomax);

for nBS = 1:96
    fit_path     = fullfile(res_path, model, 'bootstrap', num2str(nBS));
    mfs_fn   = fullfile(fit_path, 'mfs.mat');
    dps_fn   = fullfile(fit_path, 'dps.mat');
    chisq_fn   = fullfile(fit_path, 'chisq.mat');
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
%%

% dps_bsmean = dps;
% dps_bsmean = rmfield(dps_bsmean,{'mask','nii_h','s','bin'});
for nparam = 1:numel(param_map)
    eval(['dps_bsmean.' param_map{nparam} ' = mean(' param_map{nparam} ', 4);']);
    for nbin = 1:numel(dps.bin)
        eval(['dps_bsmean.bin{nbin}.' param_map{nparam} ' = mean(bin{nbin}.' param_map{nparam} ', 4);']);
    end
end

%%
figure(1), clf
axh_v = dtd_dpsbins2maps(gcf,dps_bsmean);
papersize = 3*[6 5];
set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 

fig_path = fullfile(data_path,'NII_RES',model,'maps');
eval(['print ' fig_path ' -dpdf -loose'])


%%
nbin_tot = numel(bin)+1;
for nparam = 1:numel(param_map)
    eval(['bin{nbin_tot}.' param_map{nparam} ' = ' param_map{nparam} ';']);
end
%%
bin = bin(fliplr(circshift(1:nbin_tot,[1,1])));
%%
s0max = max(dps_bsmean.s0(:));
s0thresh = 0.01;

figscale = 2;

fs = figscale*6;
lw = figscale*1;

pix_i = [9 9 9]; pix_j = [10 7 4];  pix_k = [1 1 1];

pix_col = [0         0.8    0000
         1    .3         .3
         0         0         1
         0    0.7500    0.7500];


Nbins = numel(bin);

Nparams = numel(param_map);
left = .05; bottom = .07; dleft = (1-left)/Nparams; dheight = (1-bottom)/Nbins;
width = .8*dleft; height = .9*dheight;


Nhistbins = 100;

param_xlim.s0 = s0max*[0 1];
param_xlim.mdiso = 1.2e-9*[0 1];
param_xlim.msqddelta = 1*[0 1];
param_xlim.vdiso = 2.5e-19*[0 1];
param_xlim.vsqddelta = .15*[0 1];
param_xlim.cvdisosqddelta = 1.5e-10*[-1 1];
param_xscale.s0 = s0max;
param_xscale.mdiso = 1e-9;
param_xscale.msqddelta = 1;
param_xscale.vdiso = 1e-18;
param_xscale.vsqddelta = 1;
param_xscale.cvdisosqddelta = 1e-9;

figure(2), clf

for nparam = 1:numel(param_map)
    for nbin = 1:Nbins
        dps_temp = bin{nbin};
        eval(['mapdat = dps_temp.' param_map{nparam} '/param_xscale.' param_map{nparam} ';']);
        eval(['xlim = param_xlim.' param_map{nparam} '/param_xscale.' param_map{nparam} ';']);
        axh = axes('position',[left+(nparam-1)*dleft (nbin-1)*dheight+bottom width height]);
        edges = linspace(min(xlim),max(xlim),Nhistbins);
        for npix = 1:3;
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
        if nbin > 1;
            set(axh,'XTickLabel',[])
        end
        if nparam > 1;
            set(axh,'YTickLabel',[])
        end
    end
end

papersize = figscale*[8.3 5.13];
papersize = figscale*[17.78 11];
set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 

fig_path = fullfile(data_path,'NII_RES',model,'histograms');
eval(['print ' fig_path ' -dpdf -loose'])

 