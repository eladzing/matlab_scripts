
eta=logspace(-1,log10(30),1e3);

in=4618;  % index in cata3e4 for a 1e10 stellar mass disk

fg=cata.fg(in)
beta=cata.beta(in)

fb=cata.fb(in)
xi=cata.xi(in)
bulge='hern';
Mv=cata.Mv(in)
Ms=cata.Ms(in)
qmin=cata.qmin(in)
[rvir,~,~,~]=calculate_virials('mv',cata.Mv(in));
mvms=cata.Mv(in)./cata.Ms(in);
rs=cata.rd(in)./(rvir.*1e3);
cv=cata.cv(in)
halo='nfw';



gDiskStar=eta.*bfunc(eta,1);
gDiskGas=eta.*fg.*beta^3.*bfunc(eta,beta);
gDisk=expdisk_accel(eta,'fg',fg,'beta',beta);

gBulge=bulge_accel(eta,bulge,'fb',fb,'xi',xi);
gHalo=halo_accel(eta,mvms,rs,halo,'cv',cv);

gTot=gDisk+gBulge+gHalo;


fDiskStar=gDiskStar.*expdisk_density(eta,'gas','beta',beta,'fg',1.0);
fDiskGas=gDiskGas.*expdisk_density(eta,'gas','beta',beta,'fg',1.0);
fDisk=gDisk.*expdisk_density(eta,'gas','beta',beta,'fg',1.0);

fBulge=gBulge.*expdisk_density(eta,'gas','beta',beta,'fg',1.0);
fHalo=gHalo.*expdisk_density(eta,'gas','beta',beta,'fg',1.0);
fTot=fDisk+fBulge+fHalo;

figure 
h=[];

h(2)=loglog(eta,gDiskStar,'--b','linewidth',2,...
    'DisplayName','Stellar Disk');
hold on
h(3)=loglog(eta,gDiskGas,':b','linewidth',2,...
    'DisplayName','Gas Disk');

h(1)=loglog(eta,gDisk,'-b','linewidth',2,...
    'DisplayName','Total Disk');

h(4)=loglog(eta,gBulge,'-r','linewidth',2,...
    'DisplayName','Bulge');

h(5)=loglog(eta,gHalo,'-','color',[0 0.7 0],'linewidth',2,...
    'DisplayName','Halo');

h(6)=loglog(eta,gTot,'-k','linewidth',2,...
    'DisplayName','Total');
xlim([0.1 30])


hl=legend(h);
set(hl,'Fontsize',14','Interpreter','latex','Location','SouthWest')
set(gca,'Fontsize',14,'box','on')
grid on
grid minor

xlabelmine('$r/R_d$')
ylabelmine('$g/\pi G \Sigma_s$')


% figure 
% h=[];
% 
% h(1)=loglog(eta,fDiskStar,'-b','linewidth',2,...
%     'DisplayName','Stellar Disk');
% hold on
% h(2)=loglog(eta,fDiskGas,'-r','linewidth',2,...
%     'DisplayName','Gas Disk');
% 
% h(3)=loglog(eta,fDisk,'-','color',[0 0.7 0],'linewidth',2,...
%     'DisplayName','Total Disk');
% 
% h(4)=loglog(eta,fBulge,'-c','linewidth',2,...
%     'DisplayName','Bulge');
% 
% h(5)=loglog(eta,fHalo,'-m','linewidth',2,...
%     'DisplayName','Halo');
% 
% h(6)=loglog(eta,fTot,'-k','linewidth',2,...
%     'DisplayName','Total');
% 
% 
% 
% hl=legend(h);
% set(hl,'Fontsize',14','Interpreter','latex','Location','SouthWest')
% set(gca,'Fontsize',14,'box','on')
% grid on
% grid minor
% 
% xlabelmine('$r/R_d$')
% ylabelmine('$\tilde{F}$')
% 
% %% as a function of mf
% mf=exp_disk_mass(eta,beta(1));
% 
% figure 
% h=[];
% 
% h(1)=semilogy(mf,fDiskStar,'--b','linewidth',2,...
%     'DisplayName','Stellar Disk');
% hold on
% h(2)=semilogy(mf,fDiskGas,':b','linewidth',2,...
%     'DisplayName','Gas Disk');
% 
% h(3)=semilogy(mf,fDisk,'-b','linewidth',2,...
%     'DisplayName','Total Disk');
% 
% h(4)=semilogy(mf,fBulge,'-r','linewidth',2,...
%     'DisplayName','Bulge');
% 
% h(5)=semilogy(mf,fHalo,'-','color',[0 0.7 0],'linewidth',2,...
%     'DisplayName','Halo');
% 
% h(6)=semilogy(mf,fTot,'-k','linewidth',2,...
%     'DisplayName','Total');
% 
% 
% 
% hl=legend(h);
% set(hl,'Fontsize',14','Interpreter','latex','Location','SouthWest')
% set(gca,'Fontsize',14,'box','on')
% grid on
% grid minor
% 
% xlabelmine('$M(<r)/M_{gas}$')
% ylabelmine('$\tilde{F}$')
% 
% 
