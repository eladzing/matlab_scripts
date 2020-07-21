%% Plot results for gas properties in Centrals in TNG
%% load data


bp=illustris.set_env('300');
snap=99;
global DEFAULT_MATFILE_DIR

if readFlag
    load(sprintf('%s/cooling_times_z0_TNG300.mat',DEFAULT_MATFILE_DIR))
    load(sprintf('%s/BH_energyInjection_z0_TNG300.mat',DEFAULT_MATFILE_DIR))
    
    fofs=illustris.groupcat.loadHalos(bp,snap);
    subs=illustris.groupcat.loadSubhalos(bp,snap);
    
    fofs=illustris.utils.addTvirFofs(fofs);
    
    
    subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
    units; % load general unit structure in cgs.
    
end




% identify centrals
centralMask= subsInfo.isCentral(tCoolStruct.galMask);

% get usefel stuff
galMass=tCoolStruct.galMass(tCoolStruct.galMask)';  % galaxy stellar mass
gasMass=tCoolStruct.inGal.mass;
sfr=subs.SubhaloSFRinRad(tCoolStruct.galMask)';  % sfr in galaxy
ssfr=sfr./galMass+1e-16;
sfre(:,1)=sfr./gasMass(:,1);
sfre(:,2)=sfr./gasMass(:,2);

tc=tCoolStruct.inGal.meanTcMW(:,1)';
gasTemp=tCoolStruct.inGal.meanTempMW(:,1);
gasEnt=tCoolStruct.inGal.meanEntMW(:,1);

%% black hole stuff
bhQM=bhStruct.inGal.cumEngQM;
bhRM=bhStruct.inGal.cumEngRM;

% get tvir for host fofs
tvir=fofs.Group_T_Mean200(subsInfo.hostFof+1);
tvir=tvir(tCoolStruct.galMask);
tvirMean=mk_meanMedian_bin(log10(galMass(centralMask)),log10(tvir(centralMask)),'nb',20);



%% plotting stuff
mstarLab='$\log M_\star\,[\mathrm{M_\odot}]$';
tcLab='$\log t_\mathrm{cool}\,[\mathrm{Gyr}]$';
ssfrLab='$\log \mathrm{sSFR}\,[\mathrm{yr^{-1}}] $';
sfreLab='$\log \mathrm{SFRe}\,[\mathrm{yr^{-1}}] $';
mgasLab='$\log M_\mathrm{gas}\,[\mathrm{M_\odot}]$';
tempLab='$\log T_{\mathrm{gal}} \,[\mathrm{K}]$';
entLab='$\log S_{\mathrm{gal}} \,[\mathrm{KeV\,cm^2}]$';
fgsLab='$\log M_\mathrm{gas}/M_\mathrm{star}$';
QMLab='$\log \dot{E}_\mathrm{QM}$';
RMLab='$\log \dot{E}_\mathrm{KM}$';
bhLab='$\log \dot{E}_\mathrm{KM}/\dot{E}_\mathrm{QM}$';
%cmap=brewermap(256,'RdYlBu');
cmap=brewermap(256,'*Spectral');


global simDisplayName
%% Mass vs Temp
yl=[4 8];
xl=[9 12.5];
filt=fspecial('disk',8);
xdata=log10(galMass(centralMask));
ydata=log10(gasTemp(centralMask));
popCont=plot_population_contour(xdata,ydata,'smooth',filt,'noplot','xxlim',xl,'yylim',yl);

% ssfr
figure
ax1=axes;
scatter(xdata,ydata,18,log10(ssfr(centralMask)));
hold on
errorbar(tvirMean.xMedian,tvirMean.yMean,tvirMean.yStanDev/2,'-k','linewidth',1)
hb=colorbar;barTitle(hb,ssfrLab)
colormap(cmap)
xlim(xl);ylim(yl);
grid
xlabelmine(mstarLab);
ylabelmine(tempLab);

ax2=axes;
contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','of','LineColor',[0 0 0],...
    'LevelList',10:15:100,'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];
fname=sprintf('massTemp_%s_centrals_%s','ssfr',simDisplayName);
printout_fig(gcf,fname);

% sfre
figure
ax1=axes;
scatter(xdata,ydata,18,log10(sfre(centralMask)));
hold on
errorbar(tvirMean.xMedian,tvirMean.yMean,tvirMean.yStanDev/2,'-k','linewidth',1)
hb=colorbar;barTitle(hb,sfreLab)
colormap(cmap)
xlim(xl);ylim(yl);
grid
xlabelmine(mstarLab);
ylabelmine(tempLab);

ax2=axes;
contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','of','LineColor',[0 0 0],...
    'LevelList',10:15:100,'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];

fname=sprintf('massTemp_%s_centrals_%s','sfre',simDisplayName);
printout_fig(gcf,fname);



% tcool
figure
ax1=axes;
scatter(xdata,ydata,18,log10(tc(centralMask)));
hold on
errorbar(tvirMean.xMedian,tvirMean.yMean,tvirMean.yStanDev/2,'-k','linewidth',1)
hb=colorbar;barTitle(hb,tcLab)
colormap(cmap)
xlim(xl);ylim(yl);
grid
xlabelmine(mstarLab);
ylabelmine(tempLab);

ax2=axes;
contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','of','LineColor',[0 0 0],...
    'LevelList',10:15:100,'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];

fname=sprintf('massTemp_%s_centrals_%s','tcool',simDisplayName);
printout_fig(gcf,fname);

% gasMass
fgs=tCoolStruct.inGal.mass(:,2)./galMass;
figure
ax1=axes;
scatter(xdata,ydata,18,log10(fgs(centralMask)));
hold on
errorbar(tvirMean.xMedian,tvirMean.yMean,tvirMean.yStanDev/2,'-k','linewidth',1)
hb=colorbar;barTitle(hb,fgsLab)
colormap(cmap)
xlim(xl);ylim(yl);
grid
xlabelmine(mstarLab);
ylabelmine(tempLab);

ax2=axes;
contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','of','LineColor',[0 0 0],...
    'LevelList',10:15:100,'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];


fname=sprintf('massTemp_%s_centrals_%s','mf',simDisplayName);
printout_fig(gcf,fname);

% QM
figure
ax1=axes;
scatter(xdata,ydata,18,log10(bhQM(centralMask)));
hold on
errorbar(tvirMean.xMedian,tvirMean.yMean,tvirMean.yStanDev/2,'-k','linewidth',1)
hb=colorbar;barTitle(hb,QMLab)
colormap(cmap);caxis([5 9])
xlim(xl);ylim(yl);
grid
xlabelmine(mstarLab);
ylabelmine(tempLab);

ax2=axes;
contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','of','LineColor',[0 0 0],...
    'LevelList',10:15:100,'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];

fname=sprintf('massTemp_%s_centrals_%s','QM',simDisplayName);
printout_fig(gcf,fname);

% KM
figure
ax1=axes;
scatter(xdata,ydata,18,log10(bhRM(centralMask)));
hold on
errorbar(tvirMean.xMedian,tvirMean.yMean,tvirMean.yStanDev/2,'-k','linewidth',1)
hb=colorbar;barTitle(hb,RMLab)
colormap(cmap);caxis([5 9])
xlim(xl);ylim(yl);
grid
xlabelmine(mstarLab);
ylabelmine(tempLab);

ax2=axes;
contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','of','LineColor',[0 0 0],...
    'LevelList',10:15:100,'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];

fname=sprintf('massTemp_%s_centrals_%s','KM',simDisplayName);
printout_fig(gcf,fname);

%RM/QM
figure
ax1=axes;
scatter(xdata,ydata,18,log10(bhRM(centralMask)./bhQM(centralMask)))
hold on
errorbar(tvirMean.xMedian,tvirMean.yMean,tvirMean.yStanDev/2,'-k','linewidth',1)
hb=colorbar;barTitle(hb,bhLab)
colormap(cmap);caxis([-5 0]);
xlim(xl);ylim(yl);
grid
xlabelmine(mstarLab);
ylabelmine(tempLab);

ax2=axes;
contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','of','LineColor',[0 0 0],...
    'LevelList',10:15:100,'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];

fname=sprintf('massTemp_%s_centrals_%s','bhRat',simDisplayName);
printout_fig(gcf,fname);



%% mass vs entropy
yl=[-3 2.5];
xl=[9 12.5];
filt=fspecial('disk',8);
xdata=log10(galMass(centralMask));
ydata=log10(gasEnt(centralMask));
popCont=plot_population_contour(xdata,ydata,'smooth',filt,'noplot','xxlim',xl,'yylim',yl);


% ssfr
figure
ax1=axes;
scatter(xdata,ydata,18,log10(ssfr(centralMask)));
hold on
errorbar(tvirMean.xMedian,tvirMean.yMean,tvirMean.yStanDev/2,'-k','linewidth',1)
hb=colorbar;barTitle(hb,ssfrLab)
colormap(cmap)
xlim(xl);ylim(yl);
grid
xlabelmine(mstarLab);
ylabelmine(entLab);

ax2=axes;
contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','of','LineColor',[0 0 0],...
    'LevelList',10:15:100,'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];

fname=sprintf('massEnt_%s_centrals_%s','ssfr',simDisplayName);
printout_fig(gcf,fname);

% sfre
figure
ax1=axes;
scatter(xdata,ydata,18,log10(sfre(centralMask)));
hold on
errorbar(tvirMean.xMedian,tvirMean.yMean,tvirMean.yStanDev/2,'-k','linewidth',1)
hb=colorbar;barTitle(hb,sfreLab)
colormap(cmap)
xlim(xl);ylim(yl);
grid
xlabelmine(mstarLab);
ylabelmine(entLab);

ax2=axes;
contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','of','LineColor',[0 0 0],...
    'LevelList',10:15:100,'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];

fname=sprintf('massEnt_%s_centrals_%s','sfre',simDisplayName);
printout_fig(gcf,fname);


% tcool
figure
ax1=axes;
scatter(xdata,ydata,18,log10(tc(centralMask)));
hold on
errorbar(tvirMean.xMedian,tvirMean.yMean,tvirMean.yStanDev/2,'-k','linewidth',1)
hb=colorbar;barTitle(hb,tcLab)
colormap(cmap)
xlim(xl);ylim(yl);
grid
xlabelmine(mstarLab);
ylabelmine(entLab);

ax2=axes;
contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','of','LineColor',[0 0 0],...
    'LevelList',10:15:100,'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];

fname=sprintf('massEnt_%s_centrals_%s','tcool',simDisplayName);
printout_fig(gcf,fname);

% gasMass
fgs=tCoolStruct.inGal.mass(:,2)./galMass;
figure
ax1=axes;
scatter(xdata,ydata,18,log10(fgs(centralMask)));
hold on
errorbar(tvirMean.xMedian,tvirMean.yMean,tvirMean.yStanDev/2,'-k','linewidth',1)
hb=colorbar;barTitle(hb,fgsLab)
colormap(cmap)
xlim(xl);ylim(yl);
grid
xlabelmine(mstarLab);
ylabelmine(entLab);

ax2=axes;
contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','of','LineColor',[0 0 0],...
    'LevelList',10:15:100,'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];

fname=sprintf('massEnt_%s_centrals_%s','mf',simDisplayName);
printout_fig(gcf,fname);

% QM
figure
ax1=axes;
scatter(xdata,ydata,18,log10(bhQM(centralMask)));
hold on
errorbar(tvirMean.xMedian,tvirMean.yMean,tvirMean.yStanDev/2,'-k','linewidth',1)
hb=colorbar;barTitle(hb,QMLab)
colormap(cmap);caxis([5 9])
xlim(xl);ylim(yl);
grid
xlabelmine(mstarLab);
ylabelmine(entLab);

ax2=axes;
contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','of','LineColor',[0 0 0],...
    'LevelList',10:15:100,'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];

fname=sprintf('massEnt_%s_centrals_%s','QM',simDisplayName);
printout_fig(gcf,fname);

% KM
figure
ax1=axes;
scatter(xdata,ydata,18,log10(bhRM(centralMask)));
hold on
errorbar(tvirMean.xMedian,tvirMean.yMean,tvirMean.yStanDev/2,'-k','linewidth',1)
hb=colorbar;barTitle(hb,RMLab)
colormap(cmap);caxis([5 9])
xlim(xl);ylim(yl);
grid
xlabelmine(mstarLab);
ylabelmine(entLab);

ax2=axes;
contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','of','LineColor',[0 0 0],...
    'LevelList',10:15:100,'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];

fname=sprintf('massEnt_%s_centrals_%s','KM',simDisplayName);
printout_fig(gcf,fname);


%RM/QM
figure
ax1=axes;
scatter(xdata,ydata,18,log10(bhRM(centralMask)./bhQM(centralMask)))
hold on
errorbar(tvirMean.xMedian,tvirMean.yMean,tvirMean.yStanDev/2,'-k','linewidth',1)
hb=colorbar;barTitle(hb,bhLab)
colormap(cmap);caxis([-5 0]);
xlim(xl);ylim(yl);
grid
xlabelmine(mstarLab);
ylabelmine(entLab);

ax2=axes;
contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','of','LineColor',[0 0 0],...
    'LevelList',10:15:100,'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];

fname=sprintf('massEnt_%s_centrals_%s','bhRat',simDisplayName);
printout_fig(gcf,fname);



%% add CGM stuff




