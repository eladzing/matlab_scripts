%% exact disk rps model - exploring beta & alfa

Mc=1e15;
sigma=1e7./sigma_factor('50-90');

r=0.01:0.01:60; % in units of the scale radius

%% default vales
%alfa=1; % fudge factor
fc=0.1;   % cluster gas fraction
fd=0.1;  % disk gas fraction
%Mc=1e15;  % cluster mass
delv = 337; % delta_vir
zred=0;

%% generate values of beta and alpha
nn=1000;

bt_lim=log10([0.1 10]);
%bt_lim=[0.02 5];
alf_lim=[-2 0];

lbt=bt_lim(1):diff(bt_lim)/(nn-1):bt_lim(2);
beta=10.^lbt;
lalfa=alf_lim(1):diff(alf_lim)/(nn-1):alf_lim(2);
alfa=10.^lalfa;
%% calculate pval coeeficients except for alpha dependence
etap=1;  % position in cluster in units of Rv
pval0=rps_factor_expdisk('alpha',1,'fc',fc,'fd',fd,'Mc',Mc,...
    'sigma',sigma,'delv',delv,'zred',zred,'etap',etap);

%% generate pval matrix
[~, alf]=meshgrid(beta,alfa);

pval=pval0.*alf;
clear bta alf

mstrip=zeros(nn);
rpos90=zeros(nn);
for j=1:length(lbt)
    
    mf1=exp_disk_mass(r,beta(j));
    
    fg1=disk_force_reduced(r,'beta',beta(j),'fg',fd);
    
    [fgmax,ind1]=max(fg1);
    fg=fg1(ind1:end);
    mf=mf1(ind1:end);
    
    fg_mf01=interp1(mf(mf<0.5),fg(mf<0.5),0.1); % value of fg for mf=0.1
    
    mask=(pval(:,j)>fgmax);
    mstrip(:,j)=interp1(fg,mf,pval(:,j));
    
    mstrip(mask,j)=0;
    
    %calculate map of position in cluster for 90% stripping
    rpos90(:,j)=log10(1./sqrt(fg_mf01./pval(:,j)));
    
    
end


% mass rfaction map
figure
imagesc(lbt,lalfa,100.*(1-mstrip))
set(gca,'Ydir','normal');
bar=colorbar;
set(get(bar,'Title'),'String','$\%M_{strip}$','Fontsize',12,'Interpreter','latex');
caxis([0 100])

ylabelmine('$\log(\alpha)$')
xlabelmine('$\log(\beta)$')

printout_fig(gcf,'rpsdisk_mstrip_alfabeta','v')


%
%rpos map

figure
imagesc(lbt,lalfa,rpos90)
set(gca,'Ydir','normal');
bar1=colorbar;
set(get(bar1,'Title'),'String','$r_{pos}\,[R_c]$','Fontsize',12,'Interpreter','latex');
caxis([-2 1])

ylabelmine('$\log(\alpha)$')
xlabelmine('$\log(\beta)$')

printout_fig(gcf,'rpsdisk_rpos_alfabeta','v')





