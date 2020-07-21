% show where the galaxies in the subbox are found inthe ssfr-stellar mass
% plane & in the entropy stellar mass plane. 

%% perlimiaries and data loading


sim=100;
snap=99;
bp=illustris.set_env(sim);

massThresh=10^9;

global illUnits
global simDisplayName
global DEFAULT_MATFILE_DIR

illustris.utils.set_illUnits(snap);
% load subs 

load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/fofs_subs_TNG100_z0.mat')
% fofs=illustris.groupcat.loadHalos(bp,99);
%     subs=illustris.groupcat.loadSubhalos(bp,99);
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);

load([DEFAULT_MATFILE_DIR '/cooling_times_z0_' simDisplayName '.mat'])

galEnt=tCoolStruct.inGal.meanEntMW(1,:);
cgmEnt=tCoolStruct.inCGM.meanEntMW(1,:);

ssfr=illustris.utils.calc_ssfr(subs);
smass=subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit;

mask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'centrals');

m1=mask & galEnt>0;
m2=mask & cgmEnt>0;
%% prepare contours 
filt=fspecial('disk',8);
popContS=plot_population_contour(log10(smass(mask)),log10(ssfr(mask)),'smooth',filt,'noplot');

popContG=plot_population_contour(log10(smass(m1)),log10(galEnt(m1)),'smooth',filt,'noplot');
popContC=plot_population_contour(log10(smass(m2)),log10(cgmEnt(m2)),'smooth',filt,'noplot');

%% find subbox galaxies 
bp=illustris.set_env(sim,'sub0');

global subBox 

sbmin=subBox.center'-subBox.boxSize/2;
sbmax=subBox.center'+subBox.boxSize/2;
sbMask=subs.SubhaloPos(1,:)>sbmin(1) & subs.SubhaloPos(1,:)<sbmax(1) & ...
subs.SubhaloPos(2,:)>sbmin(2) & subs.SubhaloPos(2,:)<sbmax(2) & ...
subs.SubhaloPos(3,:)>sbmin(3) & subs.SubhaloPos(3,:)<sbmax(3);


sb0Mask=sbMask & mask;

bp=illustris.set_env(sim,'sub1');

global subBox

sbmin=subBox.center'-subBox.boxSize/2;
sbmax=subBox.center'+subBox.boxSize/2;
sbMask=subs.SubhaloPos(1,:)>sbmin(1) & subs.SubhaloPos(1,:)<sbmax(1) & ...
subs.SubhaloPos(2,:)>sbmin(2) & subs.SubhaloPos(2,:)<sbmax(2) & ...
subs.SubhaloPos(3,:)>sbmin(3) & subs.SubhaloPos(3,:)<sbmax(3);


sb1Mask=sbMask & mask;

%% plot ssfr 
xl=[9 12.5];
yl=[-16.5 -8];
mstarLab='$\log M_\star\,[\mathrm{M_\odot}]$';
ssfrLab='$\log \mathrm{sSFR}\,[\mathrm{yr^{-1}}] $';
% if ~exist('cax','var')
%     cax=[min(cc) max(cc)];
% end

hf=figure;
%set(gcf,'position',[1432 421 1000 750])

ax1=axes;

plot(log10(smass(sb0Mask)),log10(ssfr(sb0Mask)),'bo','markersize',8,'linewidth',2)
hold on
plot(log10(smass(sb1Mask)),log10(ssfr(sb1Mask)),'rd','markersize',8,'linewidth',2)

hl=legend('Subbox0','Subbox1');
set(hl,'Interpreter','latex','Fontsize',14);

xlim(xl);ylim(yl);
grid
xlabelmine(mstarLab,18);
ylabelmine(ssfrLab,18);
titlemine('TNG100')
set(gca,'fontsize',14)

ax2=axes;
contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','of','LineColor',[0 0 0],...
    'LevelList',[5:10:95 100],'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];

%% plot gal Entropy

yl=[-3 2.5];

entLab=sprintf('$\\log S_{\\mathrm{%s}} \\,[\\mathrm{KeV\\,cm^2}]$','Gal');
mstarLab='$\log M_\star\,[\mathrm{M_\odot}]$';

% if ~exist('cax','var')
%     cax=[min(cc) max(cc)];
% end

hf=figure;
%set(gcf,'position',[1432 421 1000 750])

ax1=axes;

plot(log10(smass(sb0Mask)),log10(galEnt(sb0Mask)),'bo','markersize',8,'linewidth',2)
hold on
plot(log10(smass(sb1Mask)),log10(galEnt(sb1Mask)),'rd','markersize',8,'linewidth',2)

hl=legend('Subbox0','List');
set(hl,'Interpreter','latex','Fontsize',14,'Location','SouthEast');

xlim(xl);ylim(yl);
grid
xlabelmine(mstarLab,18);
ylabelmine(entLab,18);
titlemine('TNG100')
set(gca,'fontsize',14)

ax2=axes;
contour(popContG.xx,popContG.yy,popContG.popContour,'ShowText','of','LineColor',[0 0 0],...
    'LevelList',[5:10:95 100],'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];


%% plot cgm Entropy

yl=[-3 2.5];

entLab=sprintf('$\\log S_{\\mathrm{%s}} \\,[\\mathrm{KeV\\,cm^2}]$','CGM');
mstarLab='$\log M_\star\,[\mathrm{M_\odot}]$';

% if ~exist('cax','var')
%     cax=[min(cc) max(cc)];
% end

hf=figure;
%set(gcf,'position',[1432 421 1000 750])

ax1=axes;

plot(log10(smass(sb0Mask)),log10(cgmEnt(sb0Mask)),'bo','markersize',8,'linewidth',2)
hold on
plot(log10(smass(sb1Mask)),log10(cgmEnt(sb1Mask)),'rd','markersize',8,'linewidth',2)

hl=legend('Subbox0','Subbox1');
set(hl,'Interpreter','latex','Fontsize',14,'Location','SouthEast');

xlim(xl);ylim(yl);
grid
xlabelmine(mstarLab,18);
ylabelmine(entLab,18);
titlemine('TNG100')
set(gca,'fontsize',14)

ax2=axes;
contour(popContC.xx,popContC.yy,popContC.popContour,'ShowText','of','LineColor',[0 0 0],...
    'LevelList',[5:10:95 100],'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];




