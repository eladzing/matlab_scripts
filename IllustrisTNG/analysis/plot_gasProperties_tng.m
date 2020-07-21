%% Plot results for gas properties in Centrals in TNG
%% load data
global simDisplayName
if readFlag
    
    global DEFAULT_MATFILE_DIR
    load([DEFAULT_MATFILE_DIR '/cooling_times_z0_' simDisplayName '.mat'])
    load([DEFAULT_MATFILE_DIR '/BH_energyInjection_z0_' simDisplayName '.mat'])
    
        
%     bhStruct.inGal=illustris.utils.read_catalog('bh_inGal','folder','bhProps');
%     tCoolStruct.inGal=illustris.utils.read_catalog('gasProps_inGal','folder','gasProperties');
    
    
%    fofs=illustris.groupcat.loadHalos(bp,99);
%    subs=illustris.groupcat.loadSubhalos(bp,99);
    
    fofs=illustris.utils.addTvirFofs(fofs);
end
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
units; % load general unit structure in cgs.
global illUnits


global popCont
%% define what we are plotting
% identify centrals
centralMask= subsInfo.isCentral(tCoolStruct.galMask);

gasField='CGM';

%global simDisplayName 

% get usefel stuff
galMass=tCoolStruct.galMass(tCoolStruct.galMask);  % galaxy stellar mass
gasMass=tCoolStruct.(['in' gasField]).gasMass(:,tCoolStruct.galMask);
sfr=subs.SubhaloSFRinRad(tCoolStruct.galMask);  % sfr in galaxy
ssfr=illustris.utils.calc_ssfr(subs);% sfr./galMass + 10^0.5*1e-17.*10.^(0.5.*rand(size(sfr)));
ssfr=ssfr(tCoolStruct.galMask);


fgs=tCoolStruct.(['in' gasField]).gasMass(2,tCoolStruct.galMask)./galMass;
fg=fgs./(1+fgs);   % gas to total 


% 
% sfre(1,:)=sfr./gasMass(1,:);
% sfre(2,:)=sfr./gasMass(1,:);



tc=tCoolStruct.(['in' gasField]).meanTcMW(1,tCoolStruct.galMask);
gasTemp=tCoolStruct.(['in' gasField]).meanTempMW(1,tCoolStruct.galMask);
gasEnt=tCoolStruct.(['in' gasField]).meanEntMW(1,tCoolStruct.galMask);
gasMach=tCoolStruct.(['in' gasField]).medianMach(1,tCoolStruct.galMask);
gasMach2=tCoolStruct.(['in' gasField]).meanMachEW(1,tCoolStruct.galMask);
gasEDisp=tCoolStruct.(['in' gasField]).EnergyDissipation(:,tCoolStruct.galMask);
gasVDisp=tCoolStruct.(['in' gasField]).velDisp(:,tCoolStruct.galMask);
gasDens=tCoolStruct.(['in' gasField]).meanDensN(1,tCoolStruct.galMask);

galaxyMask=centralMask & gasTemp>0;
galType='centrals';


nameTag=[gasField '_' galType '_' simDisplayName];


%% black hole stuff
bhQM=bhStruct.inGal.cumEngQM(tCoolStruct.galMask);
bhRM=bhStruct.inGal.cumEngRM(tCoolStruct.galMask);
bhMass=bhStruct.inGal.bhMassMax(tCoolStruct.galMask);
%subs.SubhaloMassInRadType(illustris.partTypeNum('bh')+1,tCoolStruct.galMask).*illUnits.massUnit;


%% get tvir for host fofs
tvir=fofs.Group_T_Crit200(subsInfo.hostFof+1);
tvir=tvir(tCoolStruct.galMask);
global tvirMean
tvirMean=mk_meanMedian_bin(log10(galMass(galaxyMask)),log10(tvir(galaxyMask)),'nb',20);



%% plotting stuff
mstarLab='$\log M_\star\,[\mathrm{M_\odot}]$';
tcLab=sprintf('$\\log t_\\mathrm{cool,%s}\\,[\\mathrm{Gyr}]$',gasField);
ssfrLab='$\log \mathrm{sSFR}_{\mathrm{Gal}},[\mathrm{yr^{-1}}] $';
%sfreLab=sprintf('$\\log \\mathrm{SFRe}\\,[\mathrm{yr^{-1}}] $';
mgasLab=sprintf('$\\log M_\\mathrm{gas,%s}\\,[\\mathrm{M_\\odot}]$',gasField);
tempLab=sprintf('$\\log T_{\\mathrm{%s}} \\,[\\mathrm{K}]$',gasField);
temp2Lab=sprintf('$\\log T_{\\mathrm{%s}}/T_{\\mathrm{vir}} $',gasField);
entLab=sprintf('$\\log S_{\\mathrm{%s}} \\,[\\mathrm{KeV\\,cm^2}]$',gasField);
fgsLab=sprintf('$\\log M_\\mathrm{gas,%s}/M_\\mathrm{star}$',gasField);
fgLab=sprintf('$\\log M_\\mathrm{gas,%s}/M_\\mathrm{Baryon}$',gasField);
densLab=sprintf('$\\log n_\\mathrm{gas,%s}\\,[\\mathrm{M_\\odot/kpc^3}]$',gasField);
QMLab='$\log \dot{E}_\mathrm{QM}$';
RMLab='$\log \dot{E}_\mathrm{KM}$';
bhLab='$\log \dot{E}_\mathrm{KM}/\dot{E}_\mathrm{QM}$';
bhLab2='$\log M_{\mathrm{BH}}\,[\mathrm{M_\odot}]$';
edissLab=sprintf('$\\log \\dot{E}_\\mathrm{dis,%s}$',gasField);
sigmaLab=sprintf('$\\log \\sigma_{\\mathrm{%s}}\\,[\\mathrm{km/sec}]$',gasField);
machLab=sprintf('$\\log \\mathcal{M}_\\mathrm{%s}$',gasField);
machLab2=sprintf('$\\log \\mathcal{M}_\\mathrm{EW,%s}$',gasField);

%cmap=brewermap(256,'RdYlBu');
cmap=brewermap(256,'*Spectral');




%% Mass vs Temp
if tempFlag 
% create contours
filt=fspecial('disk',6);
xdata=log10(galMass(galaxyMask));
ydata=log10(gasTemp(galaxyMask));

popCont=plot_population_contour(xdata,ydata,'smooth',filt,'noplot');
yl=[4 8];
xl=[9 12.5];

% create tree
minLev=6;
splitParam=30;
tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);


% ssfr
plot_massTemp('tree',xdata,ydata,log10(ssfr(galaxyMask)),tre,ssfrLab,flipud(cmap),xl,yl,gasField,[-15 -9])
fname=sprintf('massTemp_%s_%s','ssfr',nameTag);
if  printFlag; printout_fig(gcf,fname,'v'); end 

% sfre
%plot_massTemp('tree',xdata,ydata,log10(sfre(galaxyMask)),sfreLab,flipud(cmap))

%fname=sprintf('massTemp_%s_%s','sfre',nameTag);
%if  printFlag; if  printFlag; printout_fig(gcf,fname,'v'); end end

% tcool
plot_massTemp('tree',xdata,ydata,log10(tc(galaxyMask)),tre,tcLab,cmap,xl,yl,gasField,[-2.5 1.5])

fname=sprintf('massTemp_%s_%s','tcool',nameTag);
if  printFlag; printout_fig(gcf,fname,'v'); end 

% gasMass

plot_massTemp('tree',xdata,ydata,log10(fg(galaxyMask)),tre,fgLab,cmap,xl,yl,gasField,[-3.5 0])

fname=sprintf('massTemp_%s_%s','mf',nameTag);
if  printFlag; printout_fig(gcf,fname,'v'); end

% gas density 
plot_massTemp('tree',xdata,ydata,log10(gasDens(galaxyMask)),tre,densLab,cmap,xl,yl,gasField,[-3 -0.5])
fname=sprintf('massTemp_%s_%s','dens',nameTag);
if  printFlag; printout_fig(gcf,fname,'v'); end

% QM
plot_massTemp('tree',xdata,ydata,log10(bhQM(galaxyMask)),tre,QMLab,cmap,xl,yl,gasField,[5 9])
fname=sprintf('massTemp_%s_%s','QM',nameTag);
if  printFlag; printout_fig(gcf,fname,'v'); end

% KM
plot_massTemp('tree',xdata,ydata,log10(bhRM(galaxyMask)),tre,RMLab,cmap,xl,yl,gasField,[5 9])

fname=sprintf('massTemp_%s_%s','KM',nameTag);
if  printFlag; printout_fig(gcf,fname,'v'); end

%RM/QM
plot_massTemp('tree',xdata,ydata,log10(bhRM(galaxyMask)./bhQM(galaxyMask)),tre,bhLab,cmap,xl,yl,gasField,[-5 0])

fname=sprintf('massTemp_%s_%s','bhRat',nameTag);
if  printFlag; printout_fig(gcf,fname,'v'); end

%BH Mass
plot_massTemp('tree',xdata,ydata,log10(bhMass(galaxyMask)),tre,bhLab2,cmap,xl,yl,gasField,[6 10])
fname=sprintf('massTemp_%s_%s','bhMass',nameTag);
if  printFlag; printout_fig(gcf,fname,'v'); end

% % Mach
% plot_massTemp('tree',xdata,ydata,log10(gasMach(galaxyMask)),tre,machLab,cmap,gasField)
% 
% fname=sprintf('massTemp_%s_%s','Mach',nameTag);
% if  printFlag; printout_fig(gcf,fname,'v'); end
% 
% % Mach ew 
% plot_massTemp('tree',xdata,ydata,log10(gasMach2(galaxyMask)),tre,machLab2,cmap,gasField)
% 
% fname=sprintf('massTemp_%s_%s','MachEW',nameTag);
% if  printFlag; printout_fig(gcf,fname,'v'); end
% 
% % e dissipation
% plot_massTemp('tree',xdata,ydata,log10(gasEDisp(2,galaxyMask)),tre,edissLab,cmap,xl,yl,gasField,[0 6])

% fname=sprintf('massTemp_%s_%s','eDisp',nameTag);
% if  printFlag; printout_fig(gcf,fname,'v'); end
% 
% % sigma - vel dispersion
% plot_massTemp('tree',xdata,ydata,log10(gasVDisp(1,galaxyMask)),tre,sigmaLab,cmap,xl,yl,gasField,[1 2.5])
% 
% fname=sprintf('massTemp_%s_%s','sigma',nameTag);
% if  printFlag; printout_fig(gcf,fname,'v'); end

end


%% mass vs entropy

if entFlag

% create contours
filt=fspecial('disk',6);
xdata=log10(galMass(galaxyMask));
ydata=log10(gasEnt(galaxyMask));
popCont=plot_population_contour(xdata,ydata,'smooth',filt,'noplot');
yl=[-3 2.5];
xl=[9 12.5];

% create tree
minLev=6;
splitParam=30;
tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);


% ssfr
plot_massEnt('tree',xdata,ydata,log10(ssfr(galaxyMask)),tre,ssfrLab,flipud(cmap),xl,yl,gasField,[-15 -9])

fname=sprintf('massEnt_%s_%s','ssfr',nameTag);
if  printFlag; printout_fig(gcf,fname,'v','v'); end

% sfre
%plot_massEnt('tree',xdata,ydata,log10(sfre(galaxyMask)),sfreLab,flipud(cmap))
%fname=sprintf('massEnt_%s_%s','sfre',nameTag);
%if  printFlag; printout_fig(gcf,fname,'v','v'); end


% tcool
plot_massEnt('tree',xdata,ydata,log10(tc(galaxyMask)),tre,tcLab,cmap,xl,yl,gasField,[-2.5 1.5])
fname=sprintf('massEnt_%s_%s','tcool',nameTag);
if  printFlag; printout_fig(gcf,fname,'v','v'); end

% gasMass

% plot_massEnt('tree',xdata,ydata,log10(fgs(galaxyMask)),tre,fgsLab,cmap,gasField)
plot_massEnt('tree',xdata,ydata,log10(fg(galaxyMask)),tre,fgLab,cmap,xl,yl,gasField,[-3.5 0])
fname=sprintf('massEnt_%s_%s','mf',nameTag);
if  printFlag; printout_fig(gcf,fname,'v','v'); end

% gas density 
plot_massEnt('tree',xdata,ydata,log10(gasDens(galaxyMask)),tre,densLab,cmap,xl,yl,gasField,[-3 -0.5])
fname=sprintf('massEnt_%s_%s','dens',nameTag);
if  printFlag; printout_fig(gcf,fname,'v'); end


% QM
plot_massEnt('tree',xdata,ydata,log10(bhQM(galaxyMask)),tre,QMLab,cmap,xl,yl,gasField,[5 9])
fname=sprintf('massEnt_%s_%s','QM',nameTag);
if  printFlag; printout_fig(gcf,fname,'v'); end

% KM
plot_massEnt('tree',xdata,ydata,log10(bhRM(galaxyMask)),tre,RMLab,cmap,xl,yl,gasField,[5 9])
fname=sprintf('massEnt_%s_%s','KM',nameTag);
if  printFlag; printout_fig(gcf,fname,'v'); end


%RM/QM
plot_massEnt('tree',xdata,ydata,log10(bhRM(galaxyMask)./bhQM(galaxyMask)),tre,bhLab,cmap,xl,yl,gasField,[-5 0])
fname=sprintf('massEnt_%s_%s','bhRat',nameTag);
if  printFlag; printout_fig(gcf,fname,'v'); end
 
%BH Mass
plot_massEnt('tree',xdata,ydata,log10(bhMass(galaxyMask)),tre,bhLab2,cmap,xl,yl,gasField,[6 10])
fname=sprintf('massEnt_%s_%s','bhMass',nameTag);
if  printFlag; printout_fig(gcf,fname,'v'); end

% Temperature 
tf=gasTemp(galaxyMask)./tvir(galaxyMask);
plot_massEnt('tree',xdata,ydata,log10(tf),tre,temp2Lab,cmap,xl,yl,gasField,[-1.5 1.5])
fname=sprintf('massEnt_%s_%s','temp',nameTag);
if  printFlag; printout_fig(gcf,fname,'v'); end




% % Mach
% plot_massEnt('tree',xdata,ydata,log10(gasMach(galaxyMask)),tre,machLab,cmap,gasField)
% 
% fname=sprintf('massEnt_%s_%s','Mach',nameTag);
% if  printFlag; printout_fig(gcf,fname,'v'); end
% 
% % Mach ew 
% plot_massEnt('tree',xdata,ydata,log10(gasMach2(galaxyMask)),tre,machLab2,cmap,gasField)
% 
% fname=sprintf('massEnt_%s_%s','MachEW',nameTag);
% if  printFlag; printout_fig(gcf,fname,'v'); end
% 
% 
% 
% % e dissipation
% plot_massEnt('tree',xdata,ydata,log10(gasEDisp(2,galaxyMask)),tre,edissLab,cmap,xl,yl,gasField,[0 6])
% 
% fname=sprintf('massEnt_%s_%s','eDisp',nameTag);
% if  printFlag; printout_fig(gcf,fname,'v'); end
% 
% % sigma - vel dispersion
% plot_massEnt('tree',xdata,ydata,log10(gasVDisp(1,galaxyMask)),tre,sigmaLab,cmap,xl,yl,gasField,[1 2.5])
% 
% fname=sprintf('massEnt_%s_%s','sigma',nameTag);
% if  printFlag; printout_fig(gcf,fname,'v'); end

end 

%% Mass vs. SSFR

% filt=fspecial('disk',6);
% xdata=log10(galMass(galaxyMask));
% ydata=log10(ssfr(galaxyMask));
% popCont=plot_population_contour(xdata,ydata,'smooth',filt,'noplot');
% xl=[9 12.5];
% yl=[-16.5 -8];
% tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
% 
% 
% % tcool
% plot_massSsfr('tree',xdata,ydata,log10(tc(galaxyMask)),tre,tcLab,cmap,[-2.5 1.5]);
% fname=sprintf('massSsfr_%s_%s','tcool',nameTag);
% if  printFlag; printout_fig(gcf,fname,'v'); end
% 
% % gasMass
% 
% plot_massSsfr('tree',xdata,ydata,log10(fgs(galaxyMask)),tre,fgsLab,cmap,[-3 0.5]);
% 
% fname=sprintf('massSsfr_%s_%s','mf',nameTag);
% if  printFlag; printout_fig(gcf,fname,'v'); end
% 
% % QM
% plot_massSsfr('tree',xdata,ydata,log10(bhQM(galaxyMask)),tre,QMLab,cmap,[5 9]);
% fname=sprintf('massSsfr_%s_%s','QM',nameTag);
% if  printFlag; printout_fig(gcf,fname,'v'); end
% 
% % KM
% plot_massSsfr('tree',xdata,ydata,log10(bhRM(galaxyMask)),tre,RMLab,cmap,[5 9]);
% fname=sprintf('massSsfr_%s_%s','KM',nameTag);
% if  printFlag; printout_fig(gcf,fname,'v'); end
% 
% 
% %RM/QM
% plot_massSsfr('tree',xdata,ydata,log10(bhRM(galaxyMask)./bhQM(galaxyMask)),tre,bhLab,cmap,[-5 0]);
% fname=sprintf('massSsfr_%s_%s','bhRat',nameTag);
% if  printFlag; printout_fig(gcf,fname,'v'); end
% 
% 
% %BH Mass
% plot_massSsfr('tree',xdata,ydata,log10(bhMass),tre,bhLab2,cmap,[6 10])
% fname=sprintf('massSsfr_%s_%s','bhMass',nameTag);
% if  printFlag; printout_fig(gcf,fname,'v'); end
% 
% 
% % Mach
% plot_massSsfr('tree',xdata,ydata,log10(gasMach(galaxyMask)),tre,machLab,cmap);
% 
% fname=sprintf('massSsfr_%s_%s','Mach',nameTag);
% if  printFlag; printout_fig(gcf,fname,'v'); end
% 
% % Mach ew 
% plot_massSsfr('tree',xdata,ydata,log10(gasMach2(galaxyMask)),tre,machLab2,cmap)
% 
% fname=sprintf('massSsfr_%s_%s','MachEW',nameTag);
% if  printFlag; printout_fig(gcf,fname,'v'); end
% 
% 
% 
% % e dissipation
% plot_massSsfr('tree',xdata,ydata,log10(gasEDisp(2,galaxyMask)),tre,edissLab,cmap,[0 6]);
% 
% fname=sprintf('massSsfr_%s_%s','eDisp',nameTag);
% if  printFlag; printout_fig(gcf,fname,'v'); end
% 
% % sigma - vel dispersion
% plot_massSsfr('tree',xdata,ydata,log10(gasVDisp(1,galaxyMask)),tre,sigmaLab,cmap,[1 2.5]);
% 
% fname=sprintf('massSsfr_%s_%s','sigma',nameTag);
% if  printFlag; printout_fig(gcf,fname,'v'); end
% 
% 
% % temp
% plot_massSsfr('tree',xdata,ydata,log10(gasTemp(galaxyMask)),tre,tempLab,cmap,[4 8]);
% 
% fname=sprintf('massSsfr_%s_%s','Temp',nameTag);
% if  printFlag; printout_fig(gcf,fname,'v'); end
% 
% 
% % Ent
% plot_massSsfr('tree',xdata,ydata,log10(gasEnt(galaxyMask)),tre,entLab,cmap,[-3 2.5]);
% 
% fname=sprintf('massSsfr_%s_%s','Ent',nameTag);
% if  printFlag; printout_fig(gcf,fname,'v'); end
% 
% 
% 
% 
% 
