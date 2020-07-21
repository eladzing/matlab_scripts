%-- 18/01/19 10:30:34 am CET --%
bp=illustris.set_env(100);
load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/cooling_times_z0_TNG100.mat')
load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/fofs_subs_TNG100_z0.mat')
load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/list_history_TNG100.mat')
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
units; % load general unit structure in cgs.
global illUnits
mask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',10^9,'centrals');
%gasEnt=tCoolStruct.inGal.meanEntMW(1,:);
gasEnt=tCoolStruct.inCGM.meanEntMW(1,:);
m1=mask & gasEnt>0;
ydata=log10(gasEnt(m1));
galMass=tCoolStruct.galMass(:)';  % galaxy stellar mass
xdata=log10(galMass(m1));
filt=fspecial('disk',8);
popCont=plot_population_contour(xdata,ydata,'smooth',filt,'noplot');

mstarLab='$\log M_\star\,[\mathrm{M_\odot}]$';
%entLab='$\log S_{\mathrm{gal}} \,[\mathrm{KeV\,cm^2}]$';
entLab=sprintf('$\\log S_{\\mathrm{%s}} \\,[\\mathrm{KeV\\,cm^2}]$','CGM');


figure
yl=[-3 2.5];
xl=[9 12.5];
cmap=brewermap(8,'Set1');

ax1=axes;
for i=1:length(galHistory)
xx=log10(galHistory(i).stellarMass);
yy=log10(galHistory(i).inCGM.meanEntMW(1,:));
%yy=log10(galHistory(i).inGal.meanEntMW(1,:));

h(i)=plot(xx,yy,':','color',cmap(i,:),'linewidth',1.5,...
'DisplayName',num2str(galHistory(i).SubfindID(1)));

if i==1
hold on
end


mm=galHistory(i).zred<=zredMax(i);

h(i)=plot(xx(mm),yy(mm),'color',cmap(i,:),'linewidth',2,...
'DisplayName',num2str(galHistory(i).SubfindID(1)));


plot(xx(1),yy(1),'o','color',cmap(i,:),'markersize',8)

plot(xx(1),yy(1),'o','color',cmap(i,:),'markersize',8)

end

hl=legend(h);
set(hl,'interpreter','latex','location','SouthEast')
xlim(xl);ylim(yl);
grid
xlabelmine(mstarLab,18);
ylabelmine(entLab,18);
set(gca,'fontsize',14)

ax2=axes;
contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','of','LineColor',[0 0 0],...
'LevelList',[5:10:95 100],'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];
