%% draw kennicutt schmidt law 

sig=10.^(6:0.01:10);

a0=2.5e-4;
a1=(2.5+0.7)*1e-4;
a2=(2.5-0.7)*1e-4;
alf=1.4;
alf1=1.4+0.15;
alf2=1.4-0.15;

td0=1/a0.*sig.^(1-alf).*10^(6*alf)./1e9; 
td1=1/a1.*sig.^(1-alf1).*10^(6*alf1)./1e9;
td2=1/a2.*sig.^(1-alf2).*10^(6*alf2)./1e9;

td1=fliplr(td1);
ll=cat(2,td2,td1);
rl=cat(2,sig,fliplr(sig));

tv=(1/sqrt(G*deltavir(0)*rho_mean(0).*Ms/Mpc^3*4*pi/3))./(yr*1e9);


figure
h=[];
patch(log10(rl),log10(ll),[0 0.65 0],'faceAlpha',0.2,'EdgeColor','none');
hold on
h(1)=loglog((sig),(td0),'-','color',[0 0.65 0],'linewidth',3,'DisplayName','Kennicutt 1998');
hold on
h(2)=loglog(xlim,(mean(dd(:,1))).*[1 1],'--b','linewidth',2,...
    'DisplayName','$\overline{t}_{\mathrm{travel}}(z=0)$');

h(3)=loglog(xlim,(mean(ddd(:,1))).*[1 1],'--r','linewidth',2,...
    'DisplayName','$\overline{t}_{\mathrm{travel}}(z=0.6)$');

h(4)=loglog(xlim,(tv).*[1 1],'-.k','linewidth',2,...
    'DisplayName','$(R_{\mathrm{vir}}/V_{\mathrm{vir}})_{z=0}$');

hl=legend(h);
set(hl,'Fontsize',14','Interpreter','latex','Location','SouthWest','box','on')
set(gca,'fontsize',14,'box','on');
grid
xlabelmine('$\log(\overline{\Sigma}_{gas})\,[\mathrm{M_{\odot}\,kpc^{-2}}]$');
ylabelmine('$\log(t_\mathrm{depl})\,[\mathrm{Gyr}]$');


