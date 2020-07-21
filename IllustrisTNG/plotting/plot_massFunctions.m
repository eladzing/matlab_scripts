%% script to plot mass function of fofs for tng 100, tng300; 
global illUnits
units;

snap=99; %z=0;
%% tng 100
% bp=illustris.set_env('100');
% 
% fofs=illustris.groupcat.loadHalos(bp,snap);
% 
% m200=fofs.Group_M_Mean200.*illUnits.massUnit;
% 
% 
% [ys,xs]=hist(log10(m200),9.5:15.5); 
% 
% x1=xs-0.5;
% x2=xs+0.5;
% x100(1:2:13)=x1;
% x100(2:2:14)=x2;
% 
% 
% %ys=fliplr(log10(cumsum(fliplr(ys))));
% ys=log10(ys);
% ys(isinf(ys))=0;
% y100(1:2:13)=ys;
% y100(2:2:14)=ys;

x100=[9 10 10 11 11 12 12 13 13 14 14 15 15 16];
y100=[6.7899 6.7899 5.0347 5.0347 4.2028 4.2028 3.3446 3.3446 2.3927 2.3927 1.3802 1.3802 0 0];
%% tng 300
bp=illustris.set_env('300');

fofs=illustris.groupcat.loadHalos(bp,snap);

m200=fofs.Group_M_Mean200.*illUnits.massUnit;


[ys,xs]=hist(log10(m200),9.5:15.5); 

x1=xs-0.5;
x2=xs+0.5;
x300(1:2:13)=x1;
x300(2:2:14)=x2;

ys=log10(ys);ys(isinf(ys))=0;
y300(1:2:13)=ys;
y300(2:2:14)=ys;


x300=[9 10 10 11 11 12 12 13 13 14 14 15 15 16];
y300=[7.1805 7.1805 6.3249 6.3249 5.4888 5.4888 4.6642 4.6642 3.6946 3.6946 2.6571 2.6571 0.7782 0.7782];


%% plot

figure
h=[];
h(1)=area(x300,y300,'FaceColor','b','linewidth',2,'FaceAlpha',0.5,'DisplayName','TNG300');
hold on
h(2)=area(x100-0.05,y100,'FaceColor','r','linewidth',2,'FaceAlpha',0.5,'DisplayName','TNG100');

set(gca,'Fontsize',14)
grid
hl=legend(h);
set(hl,'Fontsize',14,'Interpreter','latex','box','on')
xlim([9 16])
ylim([0 7.5])
xlabelmine('$\log\, M_{\mathrm{200}}\,[\mathrm{M_\odot}]$',14)
ylabelmine('$\log\, N$',14)

