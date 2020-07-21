%% draw kennicutt schmidt law 
load('/home/zinger/workProjects/matlab_scripts/cluster/mat_files/travel_times.mat')
for i=1:16
dd(i,:)=[clz0(i).dt2 clz0(i).dt1];
ddd(i,:)=[clz06(i).dt2 clz06(i).dt1];
end
tm1=min(dd(:,2));
tx1=max(dd(:,2));
tm6=min(ddd(:,2));
tx6=max(ddd(:,2));

units

sig=10.^(5:0.01:9);
xp=[sig(1) sig(end) sig(end) sig(1)];
yp1=[tm1 tm1 tx1 tx1];
yp6=[tm6 tm6 tx6 tx6];

alf=1.4;
a0=2.5e-4.*10.^(-6*alf);
eps=0.1;
tau=1; 


td0=(eps^(1-alf)-1).*sig.^(1-alf)./(alf-1)./(a0.*(1+tau))./1e9;




tv=(1/sqrt(Units.G*deltavir(0)*rho_mean(0).*Units.Ms/Units.Mpc^3*4*pi/3))./(Units.yr*1e9);


figure
h=[];
%patch(log10(rl),log10(ll),[0 0.65 0],'faceAlpha',0.2,'EdgeColor','none');

h(1)=loglog((sig),(td0),'-','color',[0 0.65 0],'linewidth',3,'DisplayName','Kennicutt 1998');
hold on
h(2)=loglog(xlim,(median(dd(:,1))).*[1 1],'--b','linewidth',2,...
    'DisplayName','$\overline{t}_{\mathrm{travel}}(z=0)$');
%fill(xp,yp1,'b','FaceAlpha',0.1);


h(3)=loglog(xlim,(median(ddd(:,1))).*[1 1],'--r','linewidth',2,...
    'DisplayName','$\overline{t}_{\mathrm{travel}}(z=0.6)$');
%fill(xp,yp6,'r','FaceAlpha',0.1);

h(4)=loglog(xlim,(tv).*[1 1],'--k','linewidth',2,...
    'DisplayName','$(R_{\mathrm{vir}}/V_{\mathrm{vir}})_{z=0}$');

ylim([0.5 100]);
hl=legend(h);
set(hl,'Fontsize',14','Interpreter','latex','Location','NorthEast','box','on')
set(gca,'fontsize',14,'box','on','Xtick',10.^(5:10));
grid
xlabelmine('$\overline{\Sigma}_{gas}\,[\mathrm{M_{\odot}\,kpc^{-2}}]$');
ylabelmine('$t_\mathrm{depl}\,[\mathrm{Gyr}]$');


