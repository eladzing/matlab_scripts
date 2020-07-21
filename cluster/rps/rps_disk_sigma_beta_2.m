%% exact disk rps model
% exploring sigma & bets for a given cluster

beta=0.1:0.05:5;

r=0.01:0.01:20; % in units of the scale radius



alfa=1; % fudge factor
fc=0.12;   % cluster gas fraction
fd=0.1;  % disk gas fraction
Mc=1e14;  % cluster mass
delv = 337; % delta_vir
zred=0;
re_fac=1.678347;

lm=14:0.1:15;
mm=10.^lm;
ls=6:0.5:12;
se=10.^ls;
%rd=1:0.1:6;

%creat sigma_e

mstrip=zeros(length(beta),length(ls),length(lm));
%pval=zeros(length(beta),length(ls),length(lm));
%sigma=zeros(length(ls),length(lm));
etap=(1).^-2;






for i=1:length(beta)
    fg1=disk_force_reducede(r,beta(i));
    mf1=exp_disk_mass(r,beta(i));
    [fgmax,ind1]=max(fg1);
    fg=fg1(ind1:end);
    mf=mf1(ind1:end);
    for j=1:length(se)
        for k=1:length(lm)
            %rd=sqrt(mm(j)/(2*pi*se(i)));
            
            pval=etap.*rps_factor('alpha',alfa,'fc',fc,'fd',fd,'beta',beta(i),...
                'Mc',mm(k),'sigma',se(j),'type','sigma_eff');
            
            %sigmae(i,j)=sigma_effective(mm(i),rd(j));
            if (pval<=fgmax)
                mstrip(i,j,k)=interp1(fg,mf,pval);
            else
                mstrip(i,j,k)=0;
            end
            
        end
    end
end



%figure
mstrip=100.*(1-mstrip);

bb=ones(size(mstrip));
ss=ones(size(mstrip));
cc=ones(size(mstrip));

for i=1:length(beta)
    bb(i,:,:)=beta(i);
end
for i=1:length(ls)
    ss(:,i,:)=ls(i);
end
for i=1:length(lm)
    cc(:,:,i)=lm(i);
end

% imagesc(ls,beta,)
% set(gca,'Ydir','normal');
% bar=colorbar;
% set(get(bar,'Title'),'String','$\%M_{strip}$','Fontsize',12,'Interpreter','latex');
%
% caxis([0 100]);
% xlabelmine('$\log(\Sigma_e),[\mathrm{M_\odot/kpc^2}]$')
% ylabelmine('$ \beta $')

% %%
% figure
% pv2=rps_factor('alpha',alfa,'fc',fc,'fd',fd,'Mc',Mc,'sigma',se,'type','sigma');
% ms=interp1(fg,mf,pv2);
% ms(pv2>fgmax)=0;
%
% semilogx(se,ms)
% xlabelmine('$\Sigma_e,[\mathrm{M_\odot/kpc^2}]$')
% ylabelmine('$\% M_{strip}$')
% % figure
% % imagesc(rd,lm,log10(sigmae))
% % set(gca,'ydir','normal');
% % bar=colorbar;
% % set(get(bar,'title'),'string','$\sigma_e$','fontsize',12,'interpreter','latex');
% %
%
% xlabelmine('$R_d\,[\mathrm{kpc}]$')
% ylabelmine('$\log(M_{disk})\,[\mathrm{M_\odot}]$')
