


mstar100=log10(M_gas1.*1e10);
mstar50=log10(M_gas.*1e10);

fg100=log10(f_gas1);
fg50=log10(f_gas);

fgm100=mk_meanMedian_bin(mstar100,fg100,'nb',8);
fgm50=mk_meanMedian_bin(mstar50,fg50,'nb',8);


cc=brewermap(8,'Set1');

figure
h(1)=plot(med.xMedian,med.yMedian,'-','color',cc(3,:),'linewidth',1.6,...
    'DisplayName','TNG100 Centrals');
hold on
plot(mstar50,fg50,'.','color',cc(1,:),'DisplayName','TNG50 Sats');


h(2)=plot(fgm50.xMedian,fgm50.yMedian,'color',cc(1,:),'linewidth',1.6,'DisplayName','TNG50 sats');

plot(mstar100,fg100,'.','color',cc(2,:),'DisplayName','TNG100 sats');
h(3)=plot(fgm100.xMedian,fgm100.yMedian,'color',cc(2,:),'linewidth',1.6,'DisplayName','TNG100 sats');

hl=legend(h);
set(hl,'Fontsize',14,'Interpreter','latex','location','SouthWest')

grid 
set(gca,'fontsize',14);

xlabelmine('log Stellar Mass $[\mathrm{M_\odot}]$');
ylabelmine('gas-stellar mass ratio');