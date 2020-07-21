%% plot cv-r=mvir relation 
figure
mv=10.^(9:0.05:15);
cv=cvir_Mvir(mv,0);
cvs=cvir_Mvir(mv,0,'random');
h=[0 0];
h(2)=semilogx(mv,cvs,'o','linewidth',1.5,'markersize',4,'DisplayName','Realization');
hold on
h(1)=semilogx(mv,cv,'-','linewidth',2,'markersize',4,'DisplayName','$c_{\mathrm{vir}}-M_{\mathrm{vir}}$ relation');

set(gca,'Fontsize',14,'box','on','Ytick',4:2:20,'Xtick',10.^(9:15))

xlim([1e9 1e15])
ylim([4 20])

xlabelmine('$\log(M_{\mathrm{vir}})\,[\mathrm{M_\odot}]$')
ylabelmine('$c_{\mathrm{vir}}$')

hl=legend(h);
set(hl,'Fontsize',14','Interpreter','latex','Location','NorthEast','box','on')