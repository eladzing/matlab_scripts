%% explore the reduced force per unit area for diferent 
% disk models 

lr=-3:0.01:log10(50);
r=10.^lr;

fr=disk_force_reduced(r,'fg',0.1,'beta',1);
frg05=disk_force_reduced(r,'fg',0.05,'beta',1);
frg2=disk_force_reduced(r,'fg',0.2,'beta',1);
frst=disk_force_reduced(r,'fg',0.0,'beta',1);

frb2=disk_force_reduced(r,'fg',0.1,'beta',2);
frg05b2=disk_force_reduced(r,'fg',0.05,'beta',2);
frg2b2=disk_force_reduced(r,'fg',0.2,'beta',2);
frstb2=disk_force_reduced(r,'fg',0.0,'beta',2);

frb05=disk_force_reduced(r,'fg',0.1,'beta',0.5);
frg05b05=disk_force_reduced(r,'fg',0.05,'beta',0.5);
frg2b05=disk_force_reduced(r,'fg',0.2,'beta',0.5);
frstb05=disk_force_reduced(r,'fg',0.0,'beta',0.5);

mf=exp_disk_mass(r,1);
mf2=exp_disk_mass(r,2);
mf05=exp_disk_mass(r,0.5);


h=[];
figure(1)
%plot
h(end+1)=loglog(r,fr,'-k','linewidth',2);
hold on
h(end+1)=loglog(r,frg05,':k','linewidth',2);
h(end+1)=loglog(r,frg2,'--k','linewidth',2);
h(end+1)=loglog(r,frst,'k-.','linewidth',2);

h(end+1)=loglog(r,frb2,'-b','linewidth',2);
h(end+1)=loglog(r,frg05b2,':b','linewidth',2);
h(end+1)=loglog(r,frg2b2,'--b','linewidth',2);
h(end+1)=loglog(r,frstb2,'b-.','linewidth',2);


h(end+1)=loglog(r,frb05,'-r','linewidth',2);
h(end+1)=loglog(r,frg05b05,':r','linewidth',2);
h(end+1)=loglog(r,frg2b05,'--r','linewidth',2);
h(end+1)=loglog(r,frstb05,'r-.','linewidth',2);


ylim([1e-3 3]);
hl=legend(h([1 2 3 4 5 9]),'$\beta=1,f_g=0.1$','$f_g=0.05$','$f_g=0.2$','stellar',...
    '$\beta=2$','$\beta=0.5$');
set(hl,'Interpreter','latex','Fontsize',11,'location','NorthEast')

set(gca,'fontsize',12,'box','on')
xlabelmine('$R/R_d$')
ylabelmine('$F_{disk}/(\pi G f_g \Sigma_s^2)$')
hold off

h=[];
figure(2)
%plot 2 
h(end+1)=semilogy(mf,fr,'-k','linewidth',2);
hold on
h(end+1)=semilogy(mf,frg05,':k','linewidth',2);
h(end+1)=semilogy(mf,frg2,'--k','linewidth',2);
h(end+1)=semilogy(mf,frst,'k-.','linewidth',2);

h(end+1)=semilogy(mf2,frb2,'-b','linewidth',2);
h(end+1)=semilogy(mf2,frg05b2,':b','linewidth',2);
h(end+1)=semilogy(mf2,frg2b2,'--b','linewidth',2);
h(end+1)=semilogy(mf2,frstb2,'b-.','linewidth',2);


h(end+1)=semilogy(mf05,frb05,'-r','linewidth',2);
h(end+1)=semilogy(mf05,frg05b05,':r','linewidth',2);
h(end+1)=semilogy(mf05,frg2b05,'--r','linewidth',2);
h(end+1)=semilogy(mf05,frstb05,'r-.','linewidth',2);


ylim([1e-3 3]);
xlim([1e-3 1]);
hl=legend(h([1 2 3 4 5 9]),'$\beta=1,f_g=0.1$','$f_g=0.05$','$f_g=0.2$','stellar',...
    '$\beta=2$','$\beta=0.5$');
set(hl,'Interpreter','latex','Fontsize',11,'location','NorthEast')

set(gca,'fontsize',12,'box','on')
xlabelmine('$M/M_{gas}$')
ylabelmine('$F_{disk}/(\pi G f_g \Sigma_s^2)$')
hold off