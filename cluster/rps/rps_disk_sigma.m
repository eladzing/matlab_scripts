%% exact disk rps model

r=0.01:0.01:20; % in units of the scale radius

mf=exp_disk_mass(r);

fg=disk_force_reducede(r);

[fgmax,ind1]=max(fg);
fg=fg(ind1:end);
mf=mf(ind1:end);

% %% perliminary figures
% figure
% plot(r,mf)
% xlabelmine('$r/r_d$');
% ylabelmine('$M/M_{sat}$');
% titlemine('Mass Profile of an Exp. Disk');
%
%
% figure
% plot(r,fg)
% xlabelmine('$r/r_d$');
% ylabelmine('unitless force');
% titlemine('force profile of an Exp. Disk');
%
% figure
% plot(mf,fg)
% xlabelmine('mass fraction');
% ylabelmine('force');
%titlemine('Mass Profile of an Exp. Disk');

%% exploring Msat and rd
alfa=1; % fudge factor
fc=0.1;   % cluster gas fraction
fd=0.1;  % disk gas fraction
Mc=1e15;  % cluster mass
delv = 337; % delta_vir
zred=0;
re_fac=1.678347;

lm=8.5:0.1:10.5;
ls=5:0.1:10;
mm=10.^lm;
se=10.^ls;
%rd=1:0.1:6;

%creat sigma_e

mstrip=zeros(length(ls),length(lm));
pval=zeros(length(ls),length(lm));
%sigma=zeros(length(ls),length(lm));
etap=(1).^-2;

for i=1:length(ls)
    for j=1:length(lm)
        rd=sqrt(mm(j)/(2*pi*se(i)));
        pval(i,j)=etap.*rps_factor('alpha',alfa,'fc',fc,'fd',fd,'Mc',Mc,...
            'Md',mm(j),'rd',rd);
        %sigmae(i,j)=sigma_effective(mm(i),rd(j));
        if (pval(i,j)<=fgmax)
            mstrip(i,j)=interp1(fg,mf,pval(i,j));
        else
            mstrip(i,j)=0;
        end
        
    end
end



figure 
imagesc(lm,ls,100.*(1-mstrip))
set(gca,'Ydir','normal');
bar=colorbar;
set(get(bar,'Title'),'String','$\%M_{strip}$','Fontsize',12,'Interpreter','latex');


ylabelmine('$\log(\Sigma_e),[\mathrm{M_\odot/kpc^2}]$')
xlabelmine('$\log(M_{disk})\,[\mathrm{M_\odot}]$')

%%
figure 
pv2=rps_factor('alpha',alfa,'fc',fc,'fd',fd,'Mc',Mc,'sigma',se','type','sigma_eff');
ms=interp1(fg,mf,pv2);
ms(pv2>fgmax)=0;

semilogx(se,ms)
xlabelmine('$\Sigma_e,[\mathrm{M_\odot/kpc^2}]$')
ylabelmine('$\% M_{strip}$')
% figure 
% imagesc(rd,lm,log10(sigmae))
% set(gca,'ydir','normal');
% bar=colorbar;
% set(get(bar,'title'),'string','$\sigma_e$','fontsize',12,'interpreter','latex');
% 
% 
% xlabelmine('$R_d\,[\mathrm{kpc}]$')
% ylabelmine('$\log(M_{disk})\,[\mathrm{M_\odot}]$')
