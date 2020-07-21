%% plot reduced force

eta=0.01:0.01:10;

beta=[1 0.5 2];
fg=[0.1 0.25];
BT=[0 0.2];
xi=3;

f1=disk_force_reduced(eta,'beta',beta(1),'fg',fg(1),'xi',xi,'BT',BT(1));
mf1=exp_disk_mass(eta,beta(1));

f2=disk_force_reduced(eta,'beta',beta(1),'fg',fg(1),'xi',xi,'BT',BT(2));
mf2=exp_disk_mass(eta,beta(1));

f3=disk_force_reduced(eta,'beta',beta(1),'fg',fg(2),'xi',xi,'BT',BT(1));
mf3=exp_disk_mass(eta,beta(1));

f4=disk_force_reduced(eta,'beta',beta(2),'fg',fg(1),'xi',xi,'BT',BT(1));
mf4=exp_disk_mass(eta,beta(2));

f5=disk_force_reduced(eta,'beta',beta(2),'fg',fg(1),'xi',xi,'BT',BT(2));
mf5=exp_disk_mass(eta,beta(2));

f6=disk_force_reduced(eta,'beta',beta(2),'fg',fg(2),'xi',xi,'BT',BT(1));
mf6=exp_disk_mass(eta,beta(2));

f7=disk_force_reduced(eta,'beta',beta(3),'fg',fg(1),'xi',xi,'BT',BT(1));
mf7=exp_disk_mass(eta,beta(3));

f8=disk_force_reduced(eta,'beta',beta(3),'fg',fg(1),'xi',xi,'BT',BT(2));
mf8=exp_disk_mass(eta,beta(3));

f9=disk_force_reduced(eta,'beta',beta(3),'fg',fg(2),'xi',xi,'BT',BT(1));
mf9=exp_disk_mass(eta,beta(3));

pt=0.15.*[1 1];

hf=figure;
lw=2.5;
h=[];
h(10)=semilogy([0 1],pt,'-.k',...
    'DisplayName','$\mathcal{P}$','linewidth',lw);
hold on
h(7)=semilogy(mf7,f7,'-',...
    'DisplayName',sprintf('$\\beta=%s,$',num2str(beta(3))),'linewidth',lw,'color',[0 0.7 0]);
h(8)=semilogy(mf8,f8,'--',...
    'DisplayName',sprintf('$\\beta=%s,$',num2str(beta(3)),num2str(BT(2)),num2str(fg(1))),'linewidth',lw,'color',[0 0.7 0]);
h(9)=semilogy(mf9,f9,':',...
    'DisplayName',sprintf('$\\beta=%s,$',num2str(beta(3)),num2str(BT(1)),num2str(fg(2))),'linewidth',lw,'color',[0 0.7 0]);
h(1)=semilogy(mf1,f1,'-b',...
    'DisplayName',sprintf('$\\beta=%s,\\,B/T=%s,\\,f_g=%s$',num2str(beta(1)),num2str(BT(1)),num2str(fg(1))),'linewidth',lw);
h(2)=semilogy(mf2,f2,'--b',...
    'DisplayName',sprintf('$\\beta=%s,\\,B/T=%s,\\,f_g=%s$',num2str(beta(1)),num2str(BT(2)),num2str(fg(1))),'linewidth',lw);
h(3)=semilogy(mf3,f3,':b',...
    'DisplayName',sprintf('$\\beta=%s,\\,B/T=%s,\\,f_g=%s$',num2str(beta(1)),num2str(BT(1)),num2str(fg(2))),'linewidth',lw);
h(4)=semilogy(mf4,f4,'-r',...
    'DisplayName',sprintf('$\\beta=%s,$',num2str(beta(2))),'linewidth',lw);
h(5)=semilogy(mf5,f5,'--r',...
    'DisplayName',sprintf('$\\beta=%s,$',num2str(beta(2))),'linewidth',lw);
h(6)=semilogy(mf6,f6,':r',...
    'DisplayName',sprintf('$\\beta=%s,$',num2str(beta(2))),'linewidth',lw);


xlabelmine('$M(\widetilde{R})/M_g$',14)
ylabelmine('$\mathcal{F}(\widetilde{R})$',14)
ylim([0.0001 10])

hl=legend(h([1 2 3 4 7 10]));
%hl=legend(h);
set(hl,'Fontsize',12,'Interpreter','latex','Location','SouthWest');
set(gca,'Fontsize',14,'box','on','Xgrid','on','ygrid','on','YMinorGrid','off')


hold off