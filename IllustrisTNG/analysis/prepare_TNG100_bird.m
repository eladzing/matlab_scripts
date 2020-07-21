%% prepare TNG1000 "trees" for background in method studies

bp=illustris.set_env(100);
snap=99; 

% load fofs and subs
load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/fofs_subs_TNG100_z0.mat')
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);

% load gas properties 
load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/gasProperties_snp99_TNG100.mat')

% units 
units; % load general unit structure in cgs.
illustris.utils.set_illUnits(snap)
global illUnits

% identify centrals
centralMask= subsInfo.isCentral(tCoolStruct.galMask);

% set parameters 
galMass=tCoolStruct.galMass(tCoolStruct.galMask);  % galaxy stellar mass
%gasTemp=tCoolStruct.(['in' gasField]).meanTempMW(1,tCoolStruct.galMask);

% sfr=subs.SubhaloSFRinRad(tCoolStruct.galMask);  % sfr in galaxy
% ssfr=illustris.utils.calc_ssfr(subs);% sfr./galMass + 10^0.5*1e-17.*10.^(0.5.*rand(size(sfr)));
% ssfr=ssfr(tCoolStruct.galMask);
% 
% ssfrThresh=1e-11;
% qMask=ssfr<=ssfrThresh;

yl=[-3 3];
xl=[9 12.5];
filt=fspecial('disk',6);    

minLev=7;
splitParam=30;

% galTree
gasEnt=tCoolStruct.inGal.meanEntMW(1,tCoolStruct.galMask);
galaxyMask=centralMask & gasEnt>0;

xdata=log10(galMass(galaxyMask));
ydata=log10(gasEnt(galaxyMask));

[galBird,~,~,~]=histogram2d(xdata,ydata,ones(size(xdata)),'len',50,'xlim',xl,'ylim',yl);
 galCont=plot_population_contour(xdata,ydata,'smooth',filt,'noplot');


%galTre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);

% cgmTree
gasEnt=tCoolStruct.inCGM.meanEntMW(1,tCoolStruct.galMask);
galaxyMask=centralMask & gasEnt>0;

xdata=log10(galMass(galaxyMask));
ydata=log10(gasEnt(galaxyMask));

[cgmBird,~,~,~]=histogram2d(xdata,ydata,ones(size(xdata)),'len',50,'xlim',xl,'ylim',yl);
cgmCont=plot_population_contour(xdata,ydata,'smooth',filt,'noplot');


%cgmTre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);

% coutTree
gasEnt=tCoolStruct.inOut.meanEntMW(1,tCoolStruct.galMask);
galaxyMask=centralMask & gasEnt>0;

xdata=log10(galMass(galaxyMask));
ydata=log10(gasEnt(galaxyMask));

[outBird,~,xxl,yyl]=histogram2d(xdata,ydata,ones(size(xdata)),'len',50,'xlim',xl,'ylim',yl);
outCont=plot_population_contour(xdata,ydata,'smooth',filt,'noplot');

%outTre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);







