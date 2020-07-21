%% exact disk rps model - exploring Mc and sigma

%% generate disk reduced force vs. mass fraction curve
beta=0.5;

r=0.01:0.01:60; % in units of the scale radius

mf1=exp_disk_mass(r,beta);

fg1=disk_force_reduced(r,'beta',beta,'fg',0.1);

[fgmax,ind1]=max(fg1);
fg=fg1(ind1:end);
mf=mf1(ind1:end);


fg_mf01=interp1(mf(mf<0.5),fg(mf<0.5),0.1); % value of fg for mf=0.1 


%% default vales 
alfa=1; % fudge factor
fc=0.1;   % cluster gas fraction
fd=0.1;  % disk gas fraction
Mc=1e15;  % cluster mass
delv = 337; % delta_vir
zred=0;

nn=200;

mc_lim=[13.5 15.5];
sig_lim=[7 11];

lm=mc_lim(1):diff(mc_lim)/(nn-1):mc_lim(2);
mm=10.^lm;
sigmac=sigma_cluster(mm);
lsg=sig_lim(1):diff(sig_lim)/(nn-1):sig_lim(2);
sg=10.^lsg;

sge=sg.*sigma_factor('half');
sg5090=sg.*sigma_factor('5090');

%% calculate pval coeeficients except for Sigma_s & Mc dependence 
etap=1;  % position in cluster in units of Rv 
pval0=rps_factor_expdisk('alpha',alfa,'fc',fc,'fd',fd,'sigma_cluster',1,...
    'sigma',1,'delv',delv,'zred',zred,'etap',etap);


%% generate sigma & sigma cluster matrix
[sgc sgs]=meshgrid(sigmac,sg);

sigrat=(sgc./sgs).^2;

clear sgc sgs

%% calculate stripping 
pval=pval0.*sigrat;

mask=(pval>fgmax);

mstrip=interp1(fg,mf,pval);
mstrip(mask)=0;

rpos90=log10(1./sqrt(fg_mf01./pval));

%% adjust observed sigma fro evolution 

%find stripping radius
rstrip=interp1(mf1,r,mstrip);


%% figures - mass rfaction maps 

% % mc vs sigam_star
% figure 
% imagesc(lm,lsg,100.*(1-mstrip))
% set(gca,'Ydir','normal','Fontsize',12);
% bar=colorbar;
% set(get(bar,'Title'),'String','$\%M_{strip}$','Fontsize',12,'Interpreter','latex');
% caxis([0 100])
% 
% ylabelmine('$\log(\Sigma_s)\,[\mathrm{M_\odot/kpc^2}]$')
% xlabelmine('$\log(M_{disk})\,[\mathrm{M_\odot}]$')

%printout_fig(gcf,'rpsdisk_mstrip_mc_sigs','v')


% mc vs sigam_eff
figure 
imagesc(lm,log10(sge),100.*(1-mstrip))
set(gca,'Ydir','normal','Fontsize',12);
bar=colorbar;
set(get(bar,'Title'),'String','$\%M_{strip}$','Fontsize',12,'Interpreter','latex');
caxis([0 100])

ylabelmine('$\log(\Sigma_{\mathrm{eff}})\,[\mathrm{M_\odot/kpc^2}]$')
xlabelmine('$\log(M_{disk})\,[\mathrm{M_\odot}]$')

%printout_fig(gcf,'rpsdisk_mstrip_mc_sige','v')

% mc vs sigam_5090
figure 
imagesc(lm,log10(sg5090),100.*(1-mstrip))
set(gca,'Ydir','normal','Fontsize',12);
bar=colorbar;
set(get(bar,'Title'),'String','$\%M_{strip}$','Fontsize',12,'Interpreter','latex');
caxis([0 100])

ylabelmine('$\log(\Sigma_{50-90})\,[\mathrm{M_\odot/kpc^2}]$')
xlabelmine('$\log(M_{disk})\,[\mathrm{M_\odot}]$')

%printout_fig(gcf,'rpsdisk_mstrip_mc_sig59','v')
%



% %% figures - rpos maps 
% 
% % Mc vs sigma_star
% 
% figure 
% imagesc(lm,lsg,rpos90)
% set(gca,'Ydir','normal');
% bar1=colorbar;
% set(get(bar1,'Title'),'String','$r_{pos}\,[R_c]$','Fontsize',12,'Interpreter','latex');
% 
% ylabelmine('$\log(\Sigma_s)\,[\mathrm{M_\odot/kpc^2}]$')
% xlabelmine('$\log(M_{disk})\,[\mathrm{M_\odot}]$')
% caxis([-2 1])
% 
% %printout_fig(gcf,'rpsdisk_rpos_mc_sigs','v')
% 
% % Mc vs sigma_eff
% figure 
% imagesc(lm,log10(sge),rpos90)
% set(gca,'Ydir','normal');
% bar1=colorbar;
% set(get(bar1,'Title'),'String','$r_{pos}\,[R_c]$','Fontsize',12,'Interpreter','latex');
% 
% ylabelmine('$\log(\Sigma_{\mathrm{eff}})\,[\mathrm{M_\odot/kpc^2}]$')
% xlabelmine('$\log(M_{disk})\,[\mathrm{M_\odot}]$')
% caxis([-2 1])
% 
% %printout_fig(gcf,'rpsdisk_rpos_mc_sige','v')
% 
% % Mc vs. Sigma_5090
% figure 
% imagesc(lm,log10(sg5090),rpos90)
% set(gca,'Ydir','normal');
% bar1=colorbar;
% set(get(bar1,'Title'),'String','$r_{pos}\,[R_c]$','Fontsize',12,'Interpreter','latex');
% 
% ylabelmine('$\log(\Sigma_{50-90})\,[\mathrm{M_\odot/kpc^2}]$')
% xlabelmine('$\log(M_{disk})\,[\mathrm{M_\odot}]$')
% caxis([-2 1])
% 
% %printout_fig(gcf,'rpsdisk_rpos_mc_sig59','v')
% 
% 
% 
% 


% figure 
% imagesc(rd,lm,log10(sigmae))
% set(gca,'Ydir','normal');
% bar=colorbar;
% set(get(bar,'Title'),'String','$\Sigma_0$','Fontsize',12,'Interpreter','latex');
% 
% 
% xlabelmine('$R_d\,[\mathrm{kpc}]$')
% ylabelmine('$\log(M_{disk})\,[\mathrm{M_\odot}]$')
