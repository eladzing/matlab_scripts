%% exact disk rps model
% we wish to find the radial position within a halo 
% for which a disk galaxy will be completely stripped 
% 

%% perliminaries
beta=2;

r=0.01:0.01:100; % in units of the scale radius
mf1=exp_disk_mass(r,beta);
fg1=disk_force_reducede(r,beta);

[fgmax,ind1]=max(fg1);
fg=fg1(ind1:end);
mf=mf1(ind1:end);

%% default values
alfa=1; % fudge factor
fc=0.1;   % cluster gas fraction
fd=0.1;  % disk gas fraction
Mc=1e15;  % cluster mass
delv = 337; % delta_vir
zred=0;
re_fac=1.678347;

% exploring sigma & Mc 
mf_crit=0.1;
fg_crit=interp1(mf(1:200),fg(1:200),mf_crit);


lm=13.5:0.05:15.5;
ls=6:0.05:10;
mm=10.^lm;
se=10.^ls; 

rstrip=zeros(length(ls),length(lm));


for i=1:length(se)
    for j=1:length(mm)
        
        pval=rps_factor_expdisk('alpha',alfa,'fc',fc,'fd',fd,'beta',beta,...
            'Mc',mm(j),'sigma',se(i),'type','sigma_eff');
        
        rstrip(i,j)=log10(sqrt(pval/fg_crit));
        
    end
end



figure 
imagesc(lm,ls,rstrip)
set(gca,'Ydir','normal');
bar=colorbar;
set(get(bar,'Title'),'String','$\log(r_p/R_c)$','Fontsize',12,'Interpreter','latex');

caxis([-2 0.5]);
ylabelmine('$\log(\Sigma_{eff}),[\mathrm{M_\odot/kpc^2}]$')
xlabelmine('$\log(M_{clust})\,[\mathrm{M_\odot}]$')
titlemine('Radius where $M_{strip}=90\%$')
hold on 
cc=[0.05 0.1 0.5 1.0 2.0 3.0];
C=contour(lm,ls,rstrip,'-k','linewidth',1.5);
clabel(C,'color','white')
hold off