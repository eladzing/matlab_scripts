bp=illustris.set_env(50);
load('C:\Users\eladz\Documents\workProjects\matlab_scripts\IllustrisTNG\matFiles\cosmic_jellyfish_objectTable.mat')
load('C:\Users\eladz\Documents\workProjects\matlab_scripts\IllustrisTNG\matFiles\jf_galProperties_CJF.mat')

myFigure

msk=objectTable.scoreWeighted>=0.8;
plot(log10(galProps.hostM200c),galProps.rpos./galProps.hostR200c,'.')
hold on
plot(log10(galProps.hostM200c(msk)),galProps.rpos(msk)./galProps.hostR200c(msk),'or')
ylim([0 1])
xlabelmine('log Halo Mass $[\mathrm{M_\odot}]$')
ylabelmine('$R/\mathrm{R_{200,c}}$')
grid
hold off
plot(log10(galProps.hostM200c(~msk)),galProps.rpos(~msk)./galProps.hostR200c(~msk),'.')
hold on
plot(log10(galProps.hostM200c(msk)),galProps.rpos(msk)./galProps.hostR200c(msk),'or')
ylim([0 1])
xlabelmine('log Halo Mass $[\mathrm{M_\odot}]$')
ylabelmine('$R/\mathrm{R_{200,c}}$')
grid
myAxis;
xlim([10.8 15])
ylim([0 1])

myFigure;
hmask=galProps.hostM200c>10^13.5;
rr=galProps.rpos./galProps.hostR200c;
histogram(rr(hmask))
msk=msk';
histogram(rr(hmask & ~msk))
hold on
histogram(rr(hmask & msk))
xlim([0 1.2])
ax=myAxis;
set(ax,'yscale','log')
legend({'Non-JF','JF'},'Interpreter','latex')
xlabelmine('$R/\mathrm{R_{200,c}}$');
ylabelmine('Number of objects')
titlemine('Radial position in Clusters $M_\mathrm{host}\ge10^{13.5}\mathrm{M_\odot}$')

hjf=histogram(rr(hmask & msk));
hjf.BinEdges
rojf=hjf.Values./hjf.BinEdges(2:end).^3;
hnjf=histogram(rr(hmask & ~msk));
ronjf=hnjf.Values./hnjf.BinEdges(2:end).^3;

myFigure;

semilogy(hnjf.BinEdges(2:end),ronjf)
xlim([0 1])
hold on
semilogy(hjf.BinEdges(2:end),rojf)
legend({'Non-JF','JF'},'Interpreter','latex')
myAxis;
xlabelmine('$R/\mathrm{R_{200,c}}$');
ylabelmine('Number Density $[R_\mathrm{200,c}^{-3}]$')

myFigure;
histogram(rr(hmask & ~msk))
hold on
histogram(rr(hmask & msk))
xlim([0 1.2])
ax=myAxis;
set(ax,'yscale','log')
legend({'Non-JF','JF'},'Interpreter','latex')
xlabelmine('$R/\mathrm{R_{200,c}}$');
ylabelmine('Number of objects')
titlemine('Radial position in Clusters $M_\mathrm{host}\ge10^{13.5}\mathrm{M_\odot}$')