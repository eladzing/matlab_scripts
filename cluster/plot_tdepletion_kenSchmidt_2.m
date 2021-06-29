%% draw kennicutt schmidt law 
load('C:\Users\eladzing\Documents\cluster\matlab\mat_files\travel_times.mat')
for i=1:16
dd(i,:)=[clz0(i).dt2 clz0(i).dt1];
ddd(i,:)=[clz06(i).dt2 clz06(i).dt1];
end


sig=10.^(5:0.01:9);


alf=1.4;
a0=2.5e-4.*10.^(-6*alf);
eps=0.1;
tau=1; 


td0=(eps^(1-alf)-1).*sig.^(1-alf)./(alf-1)./(a0.*(1+tau))./1e9;
% td0=1/a0.*sig.^(1-alf).*10^(6*alf)./1e9; 
% td1=1/a1.*sig.^(1-alf1).*10^(6*alf1)./1e9;
% td2=1/a2.*sig.^(1-alf2).*10^(6*alf2)./1e9;

% td1=fliplr(td1);
% ll=cat(2,td2,td1);
% rl=cat(2,sig,fliplr(sig));

tv=(1/sqrt(G*deltavir(0)*rho_mean(0).*Ms/Mpc^3*4*pi/3))./(yr*1e9);


figure
h=[];
%patch(log10(rl),log10(ll),[0 0.65 0],'faceAlpha',0.2,'EdgeColor','none');

h(1)=loglog((sig),(td0),'-','color',[0 0.65 0],'linewidth',3,'DisplayName','Kennicutt 1998');
hold on
h(2)=loglog(xlim,(mean(dd(:,1))).*[1 1],'--b','linewidth',2,...
    'DisplayName','$\overline{t}_{\mathrm{travel}}(z=0)$');

h(3)=loglog(xlim,(mean(ddd(:,1))).*[1 1],'--r','linewidth',2,...
    'DisplayName','$\overline{t}_{\mathrm{travel}}(z=0.6)$');

h(4)=loglog(xlim,(tv).*[1 1],'-.k','linewidth',2,...
    'DisplayName','$(R_{\mathrm{vir}}/V_{\mathrm{vir}})_{z=0}$');

ylim([0.1 25]);
hl=legend(h);
set(hl,'Fontsize',14','Interpreter','latex','Location','SouthWest','box','on')
set(gca,'fontsize',14,'box','on','Xtick',10.^(5:10));
grid
xlabelmine('$\log(\overline{\Sigma}_{gas})\,[\mathrm{M_{\odot}\,kpc^{-2}}]$');
ylabelmine('$\log(t_\mathrm{depl})\,[\mathrm{Gyr}]$');


