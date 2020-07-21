%% Plot results for gas properties in Centrals in TNG

%% load data

if readFlag
   global DEFAULT_MATFILE_DIR
    load([DEFAULT_MATFILE_DIR '/cooling_times_z0_TNG100.mat'])
    load([DEFAULT_MATFILE_DIR '/BH_energyInjection_z0_TNG100.mat'])
    
%     list={'Gal' 'CGM' 'Sub'};
%    for i=1:length(list)
%        
%        bhStruct.(['in' list{i}])=illustris.utils.read_catalog(['bh_in' list{i}],'folder','bhProps');
%        tCoolStruct.(['in' list{i}])=illustris.utils.read_catalog(['gasProps_in' list{i}],'folder','gasProperties');
%    end
    
    
    
    global DRACOFLAG
    if DRACOFLAG
        fofs=illustris.groupcat.loadHalos(bp,99);
        subs=illustris.groupcat.loadSubhalos(bp,99);
    else
        load([DEFAULT_MATFILE_DIR '/tng100_z0_fofs_subs.mat'])
    end
    fofs=illustris.utils.addTvirFofs(fofs);
end


subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
units; % load general unit structure in cgs.


%% define what we are plotting 
% identify centrals
centralMask= subsInfo.isCentral(tCoolStruct.galMask);

global simDisplayName
gasField='CGM'

% get useful stuff
galMass=tCoolStruct.galMass(tCoolStruct.galMask);  % galaxy stellar mass
sfr=subs.SubhaloSFRinRad(tCoolStruct.galMask);  % sfr in galaxy
ssfr=sfr./galMass + 10^0.5*1e-17.*10.^(0.5.*rand(size(sfr)));


gasMass=tCoolStruct.(['in' gasField]).gasMass(:,tCoolStruct.galMask);
tc=tCoolStruct.(['in' gasField]).meanTcMW(1,tCoolStruct.galMask);
gasTemp=tCoolStruct.(['in' gasField]).meanTempMW(1,tCoolStruct.galMask);
gasEnt=tCoolStruct.(['in' gasField]).meanEntMW(1,tCoolStruct.galMask);
gasMach=tCoolStruct.(['in' gasField]).medianMach(1,tCoolStruct.galMask);
gasMach2=tCoolStruct.(['in' gasField]).meanMachEW(1,tCoolStruct.galMask);
gasEDisp=tCoolStruct.(['in' gasField]).EnergyDissipation(:,tCoolStruct.galMask);
gasVDisp=tCoolStruct.(['in' gasField]).velDisp(:,tCoolStruct.galMask);
sfre(1,:)=sfr./gasMass(1,:);
sfre(2,:)=sfr./gasMass(1,:);

galaxyMask=centralMask & gasTemp>0;
galType='centrals'


nameTag=[gasField '_' galType '_' simDisplayName];


%% black hole stuff
bhQM=bhStruct.inGal.cumEngQM(tCoolStruct.galMask);
bhRM=bhStruct.inGal.cumEngRM(tCoolStruct.galMask);

%% get tvir for host fofs
tvir=fofs.Group_T_Mean200(subsInfo.hostFof+1);
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
entLab=sprintf('$\\log S_{\\mathrm{%s}} \\,[\\mathrm{KeV\\,cm^2}]$',gasField);
fgsLab=sprintf('$\\log M_\\mathrm{gas,%s}/M_\\mathrm{star}$',gasField);
QMLab='$\log \dot{E}_\mathrm{QM}$';
RMLab='$\log \dot{E}_\mathrm{KM}$';
bhLab='$\log \dot{E}_\mathrm{KM}/\dot{E}_\mathrm{QM}$';
edissLab=sprintf('$\\log \\dot{E}_\\mathrm{dis,%s}$',gasField);
sigmaLab=sprintf('$\\log \\sigma_{\\mathrm{%s}}\\,[\\mathrm{km/sec}]$',gasField);
machLab=sprintf('$\\log \\mathcal{M}_\\mathrm{%s}$',gasField);
machLab2=sprintf('$\\log \\mathcal{M}_\\mathrm{EW,%s}$',gasField);







%cmap=brewermap(256,'RdYlBu');
cmap=brewermap(256,'*Spectral');



%% Mass vs Temp

filt=fspecial('disk',8);
xdata=log10(galMass(galaxyMask));
ydata=log10(gasTemp(galaxyMask));
global popCont
popCont=plot_population_contour(xdata,ydata,'smooth',filt,'noplot');
yl=[4 8];
xl=[9 12.5];

% create tree 
minLev=6;
splitParam=50;
tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);

% ssfr
plot_massTemp('tree',xdata,ydata,log10(ssfr(galaxyMask)),tre,ssfrLab,flipud(cmap),gasField)
fname=sprintf('massTemp_%s_%s','ssfr',nameTag);
printout_fig(gcf,fname);

% sfre
%plot_massTemp('tree',xdata,ydata,log10(sfre(galaxyMask)),sfreLab,flipud(cmap),gasField)

%fname=sprintf('massTemp_%s_%s','sfre',nameTag);
%printout_fig(gcf,fname);

% tcool
plot_massTemp('tree',xdata,ydata,log10(tc(galaxyMask)),tre,tcLab,cmap,gasField,[-3 1.5])

fname=sprintf('massTemp_%s_%s','tcool',nameTag);
printout_fig(gcf,fname);

% gasMass
fgs=tCoolStruct.(['in' gasField]).gasMass(2,tCoolStruct.galMask)./galMass;
plot_massTemp('tree',xdata,ydata,log10(fgs(galaxyMask)),tre,fgsLab,cmap,gasField)

fname=sprintf('massTemp_%s_%s','mf',nameTag);
printout_fig(gcf,fname);

% QM

plot_massTemp('tree',xdata,ydata,log10(bhQM(galaxyMask)),tre,QMLab,cmap,gasField,[5 9])
fname=sprintf('massTemp_%s_%s','QM',nameTag);
printout_fig(gcf,fname);

% KM
plot_massTemp('tree',xdata,ydata,log10(bhRM(galaxyMask)),tre,RMLab,cmap,gasField,[5 9])

fname=sprintf('massTemp_%s_%s','KM',nameTag);
printout_fig(gcf,fname);

%RM/QM
plot_massTemp('tree',xdata,ydata,log10(bhRM(galaxyMask)./bhQM(galaxyMask)),tre,bhLab,cmap,gasField,[-5 0])

fname=sprintf('massTemp_%s_%s','bhRat',nameTag);
printout_fig(gcf,fname);


% Mach
plot_massTemp('tree',xdata,ydata,log10(gasMach(galaxyMask)),tre,machLab,cmap,gasField)

fname=sprintf('massTemp_%s_%s','Mach',nameTag);
printout_fig(gcf,fname);

% Mach ew 
plot_massTemp('tree',xdata,ydata,log10(gasMach2(galaxyMask)),tre,machLab2,cmap,gasField)

fname=sprintf('massTemp_%s_%s','MachEW',nameTag);
printout_fig(gcf,fname);


% e dissipation
plot_massTemp('tree',xdata,ydata,log10(gasEDisp(2,galaxyMask)),tre,edissLab,cmap,gasField,[0 6])

fname=sprintf('massTemp_%s_%s','eDisp',nameTag);
printout_fig(gcf,fname);

% sigma - vel dispersion
plot_massTemp('tree',xdata,ydata,log10(gasVDisp(1,galaxyMask)),tre,sigmaLab,cmap,gasField,[1 2.5])

fname=sprintf('massTemp_%s_%s','sigma',nameTag);
printout_fig(gcf,fname);


%% mass vs entropy
% create contours
filt=fspecial('disk',8);
xdata=log10(galMass(galaxyMask));
ydata=log10(gasEnt(galaxyMask));
popCont=plot_population_contour(xdata,ydata,'smooth',filt,'noplot');
yl=[-3 2.5];
xl=[9 12.5];

% create tree
minLev=6;
splitParam=20;
tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);


% ssfr
plot_massEnt('tree',xdata,ydata,log10(ssfr(galaxyMask)),tre,ssfrLab,flipud(cmap),gasField)

fname=sprintf('massEnt_%s_%s','ssfr',nameTag);
printout_fig(gcf,fname);

% sfre
%plot_massEnt('tree',xdata,ydata,log10(sfre(galaxyMask)),sfreLab,flipud(cmap),gasField)
%fname=sprintf('massEnt_%s_%s','sfre',nameTag);
%printout_fig(gcf,fname);


% tcool
plot_massEnt('tree',xdata,ydata,log10(tc(galaxyMask)),tre,tcLab,cmap,gasField,[-3 1.5])
fname=sprintf('massEnt_%s_%s','tcool',nameTag);
printout_fig(gcf,fname);

% gasMass

plot_massEnt('tree',xdata,ydata,log10(fgs(galaxyMask)),tre,fgsLab,cmap,gasField)

fname=sprintf('massEnt_%s_%s','mf',nameTag);
printout_fig(gcf,fname);

% QM
plot_massEnt('tree',xdata,ydata,log10(bhQM(galaxyMask)),tre,QMLab,cmap,gasField,[5 9])
fname=sprintf('massEnt_%s_%s','QM',nameTag);
printout_fig(gcf,fname);

% KM
plot_massEnt('tree',xdata,ydata,log10(bhRM(galaxyMask)),tre,RMLab,cmap,gasField,[5 9])
fname=sprintf('massEnt_%s_%s','KM',nameTag);
printout_fig(gcf,fname);


%RM/QM
plot_massEnt('tree',xdata,ydata,log10(bhRM(galaxyMask)./bhQM(galaxyMask)),tre,bhLab,cmap,gasField,[-5 0])
fname=sprintf('massEnt_%s_%s','bhRat',nameTag);
printout_fig(gcf,fname);


% Mach
plot_massEnt('tree',xdata,ydata,log10(gasMach(galaxyMask)),tre,machLab,cmap,gasField)

fname=sprintf('massEnt_%s_%s','Mach',nameTag);
printout_fig(gcf,fname);

% Mach ew 
plot_massEnt('tree',xdata,ydata,log10(gasMach2(galaxyMask)),tre,machLab2,cmap,gasField)

fname=sprintf('massEnt_%s_%s','MachEW',nameTag);
printout_fig(gcf,fname);

% e dissipation
plot_massEnt('tree',xdata,ydata,log10(gasEDisp(2,galaxyMask)),tre,edissLab,cmap,gasField,[0 6])

fname=sprintf('massEnt_%s_%s','eDisp',nameTag);
printout_fig(gcf,fname);

% sigma - vel dispersion
plot_massEnt('tree',xdata,ydata,log10(gasVDisp(1,galaxyMask)),tre,sigmaLab,cmap,gasField,[1 2.5])

fname=sprintf('massEnt_%s_%s','sigma',nameTag);
printout_fig(gcf,fname);


%% Mass vs. SSFR

filt=fspecial('disk',8);
xdata=log10(galMass(galaxyMask));
ydata=log10(ssfr(galaxyMask));
popCont=plot_population_contour(xdata,ydata,'smooth',filt,'noplot');
xl=[9 12.5];
yl=[-16.5 -8];
tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);


% tcool
plot_massSsfr('tree',xdata,ydata,log10(tc(galaxyMask)),tre,tcLab,cmap,[-3 1.5]);
fname=sprintf('massSsfr_%s_%s','tcool',nameTag);
printout_fig(gcf,fname);

% gasMass

plot_massSsfr('tree',xdata,ydata,log10(fgs(galaxyMask)),tre,fgsLab,cmap,[-3 0.5]);

fname=sprintf('massSsfr_%s_%s','mf',nameTag);
printout_fig(gcf,fname);


% Mach
plot_massSsfr('tree',xdata,ydata,log10(gasMach(galaxyMask)),tre,machLab,cmap);

fname=sprintf('massSsfr_%s_%s','Mach',nameTag);
printout_fig(gcf,fname);

% Mach ew 
plot_massSsfr('tree',xdata,ydata,log10(gasMach2(galaxyMask)),tre,machLab2,cmap,gasField)

fname=sprintf('massSsfr_%s_%s','MachEW',nameTag);
printout_fig(gcf,fname);



% e dissipation
plot_massSsfr('tree',xdata,ydata,log10(gasEDisp(2,galaxyMask)),tre,edissLab,cmap,[0 6]);

fname=sprintf('massSsfr_%s_%s','eDisp',nameTag);
printout_fig(gcf,fname);

% sigma - vel dispersion
plot_massSsfr('tree',xdata,ydata,log10(gasVDisp(1,galaxyMask)),tre,sigmaLab,cmap,[1 2.5]);

fname=sprintf('massSsfr_%s_%s','sigma',nameTag);
printout_fig(gcf,fname);


% temp
plot_massSsfr('tree',xdata,ydata,log10(gasTemp(galaxyMask)),tre,tempLab,cmap,[4 8]);

fname=sprintf('massSsfr_%s_%s','Temp',nameTag);
printout_fig(gcf,fname);


% Ent
plot_massSsfr('tree',xdata,ydata,log10(gasEnt(galaxyMask)),tre,entLab,cmap,[-3 2.5]);

fname=sprintf('massSsfr_%s_%s','Ent',nameTag);
printout_fig(gcf,fname);






