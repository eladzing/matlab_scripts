figure 
cc= brewermap(8,'Set1'); %distinguishable_colors(3);
h=[];

h(1)=semilogy(sp1.spec(:,1),sp1.spec(:,2),'Color',cc(1,:),'DisplayName',sprintf('$Z=%3.2f$',sp1.abun));
hold on
h(2)=semilogy(sp2.spec(:,1),sp2.spec(:,2),'Color',cc(2,:),'DisplayName',sprintf('$Z=%3.2f$',sp2.abun));
h(3)=semilogy(sp3.spec(:,1),sp3.spec(:,2),'Color',cc(3,:),'DisplayName',sprintf('$Z=%3.2f$',sp3.abun));
h(4)=semilogy(sp4.spec(:,1),sp4.spec(:,2),'Color',cc(4,:),'DisplayName',sprintf('$Z=%3.2f$',sp4.abun));
h(5)=semilogy(sp5.spec(:,1),sp5.spec(:,2),'Color',cc(5,:),'DisplayName',sprintf('$Z=%3.2f$',sp5.abun));
h(6)=semilogy(sp6.spec(:,1),sp6.spec(:,2),'Color',cc(6,:),'DisplayName',sprintf('$Z=%3.2f$',sp6.abun));

hl=legend(h);
set(hl,'Interpreter','latex','Fontsize',14);
xlabelmine('$\mathrm{keV}$')
ylabelmine('$\mathrm{counts\,s^{-1}\,keV^{-1}\,cm^{3}}$');

figure 
abun=[sp1.abun sp2.abun sp3.abun sp4.abun sp5.abun sp6.abun];
coounts=[sum(sp1.spec(:,2)) sum(sp2.spec(:,2)) sum(sp3.spec(:,2)) sum(sp4.spec(:,2)) sum(sp5.spec(:,2)) sum(sp6.spec(:,2))]./sum(sp4.spec(:,2));

plot(abun,coounts,'o')
xlabelmine('$Z\,[\mathrm{Z_{\odot}}]$');
ylabelmine('$\mathrm{counts\,s^{-1}\,keV^{-1}\,cm^{3}}$');

%% plot redshift dependence 

figure 

h=[];

h(1)=semilogy(spz(1).spec(:,1),spz(1).spec(:,2),'Color',cc(1,:),'DisplayName',sprintf('$z=%3.2f$',spz(1).zred_emit));
hold on
h(2)=semilogy(spz(2).spec(:,1),spz(2).spec(:,2),'Color',cc(2,:),'DisplayName',sprintf('$z=%3.2f$',spz(2).zred_emit));
h(3)=semilogy(spz(3).spec(:,1),spz(3).spec(:,2),'Color',cc(3,:),'DisplayName',sprintf('$z=%3.2f$',spz(3).zred_emit));
h(4)=semilogy(spz(4).spec(:,1),spz(4).spec(:,2),'Color',cc(4,:),'DisplayName',sprintf('$z=%3.2f$',spz(4).zred_emit));
h(5)=semilogy(spz(5).spec(:,1),spz(5).spec(:,2),'Color',cc(5,:),'DisplayName',sprintf('$z=%3.2f$',spz(5).zred_emit));
h(6)=semilogy(spz(6).spec(:,1),spz(6).spec(:,2),'Color',cc(6,:),'DisplayName',sprintf('$z=%3.2f$',spz(6).zred_emit));

hl=legend(h);
set(hl,'Interpreter','latex','Fontsize',14);
xlabelmine('$\mathrm{keV}$')
ylabelmine('$\mathrm{counts\,s^{-1}\,keV^{-1}\,cm^{3}}$');

figure 
abun=[spz(1).zred_emit spz(2).zred_emit spz(3).zred_emit spz(4).zred_emit spz(5).zred_emit spz(6).zred_emit];
coounts=[sum(spz(1).spec(:,2)) sum(spz(2).spec(:,2)) sum(spz(3).spec(:,2)) sum(spz(4).spec(:,2)) sum(spz(5).spec(:,2)) sum(spz(6).spec(:,2))]./sum(spz(3).spec(:,2));

plot(abun,coounts,'o')
xlabelmine('$z$');
ylabelmine('$\mathrm{counts\,s^{-1}\,keV^{-1}\,cm^{3}}$');