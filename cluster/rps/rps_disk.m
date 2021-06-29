%% exact disk rps model  - exploring Msat and rd

%% generate disk reduced force vs. mass fraction curve
beta=1;
fb=0.25;
fgs=0.1;

r=0.01:0.01:60; % in units of the scale radius

mf1=exp_disk_mass(r,beta);

fg1=disk_force_reduced(r,'beta',beta,'fg',fgs,'fb',fb);

[fgmax,ind1]=max(fg1);
fg=fg1(ind1:end);
mf=mf1(ind1:end);


fg_mf01=interp1(mf(mf<0.5),fg(mf<0.5),0.5); % value of fg for mf=0.1

%% default vales 
alfa=1; % fudge factor
fc=0.1;   % cluster gas fraction
fd=fgs;  % disk gas fraction
Mc=1e15;  % cluster mass
zred=0;
delv = deltavir(zred); % delta_vir


%% generate values of msat and rscale
nn=1000;

ms_lim=[9 11];
rd_lim=[1 6];

lm=ms_lim(1):diff(ms_lim)/(nn-1):ms_lim(2);
mm=10.^lm;
rd=rd_lim(1):diff(rd_lim)/(nn-1):rd_lim(2);

%% calculate pval coeeficients except for Sigma_s dependence 
etap=0.5;  % position in cluster in units of Rv 
pval0=rps_factor_expdisk('alpha',alfa,'fc',fc,'fd',fd,'Mc',Mc,...
    'sigma',1,'delv',delv,'zred',zred,'etap',etap);

%% generate sigma matrix
[mmg rdg]=meshgrid(mm,rd);
sigma=mmg./(2*pi*rdg.^2);

clear mmg rdg

sigmae=sigma.*sigma_factor('eff');

%% calculate stripping 
pval=pval0./sigma.^2;

mask=(pval>fgmax);

mstrip=interp1(fg,mf,pval);
mstrip(mask)=0;

rpos90=log10(1./sqrt(fg_mf01./pval));

%% figures
% mass rfaction map
figure
imagesc(lm,rd,100.*(1-mstrip))
set(gca,'Ydir','normal','Fontsize',12);
bar=colorbar;
set(get(bar,'Title'),'String','$\%M_{strip}$','Fontsize',12,'Interpreter','latex');


ylabelmine('$R_d\,[\mathrm{kpc}]$')
xlabelmine('$\log(M_{disk})\,[\mathrm{M_\odot}]$')

hold on
cc=6:0.2:10;
C=contour(lm,rd,log10(sigmae),cc,'-k','linewidth',1.5);
clabel(C,'color','white')
hold off
caxis([0 100])

%printout_fig(gcf,'rpsdisk_mstrip_mdrd','v')
%
%rpos map

figure
imagesc(lm,rd,rpos90)
set(gca,'Ydir','normal','Fontsize',12);
bar1=colorbar;
set(get(bar1,'Title'),'String','$r_{pos}\,[R_c]$','Fontsize',12,'Interpreter','latex');


ylabelmine('$R_d\,[\mathrm{kpc}]$')
xlabelmine('$\log(M_{disk})\,[\mathrm{M_\odot}]$')

hold on
cc=6:0.2:10;
C=contour(lm,rd,log10(sigmae),cc,'-k','linewidth',1.5);
clabel(C,'color','white')
hold off

caxis([-2 1])

%printout_fig(gcf,'rpsdisk_rpos_mdrd','v')





% figure
% imagesc(rd,lm,log10(sigmae))
% set(gca,'Ydir','normal');
% bar=colorbar;
% set(get(bar,'Title'),'String','$\Sigma_0$','Fontsize',12,'Interpreter','latex');
%
%
% xlabelmine('$R_d\,[\mathrm{kpc}]$')
% ylabelmine('$\log(M_{disk})\,[\mathrm{M_\odot}]$')
