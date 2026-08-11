%% This Script generates interesting gas properties for the subfind catalogs
% and saves them as catalogs for TNG


%% set framework
% the simulation, snapshot, and environment should be set prior to running
% the script

illustris.utils.set_illUnits(snap);



global illUnits
global DEFAULT_MATFILE_DIR
global simDisplayName
global cosmoStruct


if ~exist('readFlag','var')
    readFlag=true;
end

% read FOF and SUBFIND data, as well as free-fall time profiles
if readFlag
    fprintf(' *** Reading data *** \n');
    [subs, fofs, subsInfo]=illustris.loadFofSub(snap);

    % fofs=illustris.groupcat.loadHalos(bp,snap);
    % subs=illustris.groupcat.loadSubhalos(bp,snap);
    readFlag=false;
    fofs = illustris.utils.addTvirFofs( fofs);
    %load([DEFAULT_MATFILE_DIR '/freeFallTime_profiles_snp' num2str(illUnits.snap) '_' simDisplayName '.mat']);
end

% assign a lower stellar mass threshold based on the TNG box.
if contains(simDisplayName,'100') || contains(simDisplayName,'300')
    massThresh=10^9; % threshold for *stellar* mass
elseif contains(simDisplayName,'50')
    massThresh=10^(8.3); % threshold for *stellar* mass
    massThreshTop=10^(9.5);
else
    error('mass threshold not set  - could not identify simulation: %s \n',simDisplayName);
end
fprintf('Setting stellar mass threshold to: %0.1e solar mass \n',massThresh);




%% select galaxies
M200c=double(fofs.Group_M_Crit200(subsInfo.hostFof+1)).*illUnits.massUnit;
R200c=double(fofs.Group_R_Crit200(subsInfo.hostFof+1));  %.*illUnits.lengthUnit;
T200c=double(fofs.Group_T_Crit200(subsInfo.hostFof+1));
K200c=mass_entropy_relation(M200c,'zred',illUnits.zred,'cosmo',cosmoStruct,...
    'delta',200,'rhotype','crit');
massAllGals=illustris.utils.get_stellar_mass(subs,'gal');

% this mask selects all galaxies with dm component & stars, above *stellar* mass limit whose host has virials parameters
galMaskDwarf=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'massTop',massThreshTop,'snap',snap,'gas','centrals');
galMaskAll=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'gas','centrals');

dwarfIndx=find(galMaskDwarf);
allIndx=find(galMaskAll);

%% find isolation condition in 3D space 

isoThresh=5;

nneib3=illustris.utils.find_k_nearest_neighbor_3D(subs.SubhaloPos(:,galMaskAll),1,'qp',subs.SubhaloPos(:,galMaskDwarf));
nneib2=illustris.utils.find_k_nearest_neighbor_2D_vel(subs.SubhaloPos(:,galMaskAll),subs.SubhaloVel(:,galMaskAll),  ...
    1,300,'qp',subs.SubhaloPos(:,galMaskDwarf),subs.SubhaloVel(:,galMaskDwarf));

indx3=allIndx(nneib3.indx);
rnorm=max(R200c(dwarfIndx),R200c(indx3));
isolatedMask3=nneib3.distance./rnorm>=isoThresh;


%% find isolation condition in 2D-V space 
for i=1:3
    rnorm=max(R200c(dwarfIndx),R200c(allIndx(nneib2.indx(i,:))));
    isolatedMask2(i,:)=nneib2.distance(i,:)./rnorm>=isoThresh;
    isolatedMask22(i,:)=nneib2.distance(i,:)./rnorm>=isoThresh*sqrt(2/3);
    isolatedMask23(i,:)=nneib2.distance(i,:)./rnorm>=isoThresh*0.6;
end

%% compare isolatation criteria 
myFigure;
i=1;
cdfplot(log10(nneib2.distance(i,:)./rnorm));
hold on
i=2;
cdfplot(log10(nneib2.distance(i,:)./rnorm));
i=3;
cdfplot(log10(nneib2.distance(i,:)./rnorm));
cdfplot(log10(nneib3.distance./max(R200c(dwarfIndx),R200c(allIndx(nneib3.indx)))));

%% compare isolatation criteria 
% len=length(nneib3.distance);
% myFigure;
% i=1;
% hst1=histogram(log10(nneib2.distance(i,:)./rnorm),linspace(-1,log10(50),51),'Normalization','cumcount');
% hold on
% i=2;
% hst2=histogram(log10(nneib2.distance(i,:)./rnorm),linspace(-1,log10(50),51),'Normalization','cumcount');
% i=3;
% hst3=histogram(log10(nneib2.distance(i,:)./rnorm),linspace(-1,log10(50),51),'Normalization','cumcount');
% hst4=histogram(log10(nneib3.distance./max(R200c(indxBase),R200c(indx2(nneib3.indx)))),linspace(-1,log10(50),51),'Normalization','cumcount');
% hst41=histogram(log10(sqrt(2/3).*nneib3.distance./max(R200c(indxBase),R200c(indx2(nneib3.indx)))),linspace(-1,log10(50),51),'Normalization','cumcount');

%% compare isolatation criteria 
len=length(nneib3.distance);
myFigure;
i=1;
hst1=histogram((nneib2.distance(i,:)./rnorm),linspace(0,(50),101),'Normalization','cumcount');
hold on
i=2;
hst2=histogram((nneib2.distance(i,:)./rnorm),linspace(0,(50),101),'Normalization','cumcount');
i=3;
hst3=histogram((nneib2.distance(i,:)./rnorm),linspace(0,(50),101),'Normalization','cumcount');
hst4=histogram((nneib3.distance./max(R200c(dwarfIndx),R200c(allIndx(nneib3.indx)))),linspace(0,(50),101),'Normalization','cumcount');
hst41=histogram((sqrt(2/3).*nneib3.distance./max(R200c(dwarfIndx),R200c(allIndx(nneib3.indx)))),linspace(-1,(50),51),'Normalization','cumcount');

%%

colors=brewermap(9,'Set1');
myFigure;
h=[];
h(end+1)=plot(hst1.BinEdges(2:end),len-hst1.Values,'color',colors(1,:),'LineWidth',1.5,...
    'DisplayName','2Dx');
hold on
h(end+1)=plot(hst2.BinEdges(2:end),len-hst2.Values,'color',colors(2,:),'LineWidth',1.5,...
    'DisplayName','2Dy');
h(end+1)=plot(hst3.BinEdges(2:end),len-hst3.Values,'color',colors(3,:),'LineWidth',1.5,...
    'DisplayName','2Dz');
h(end+1)=plot(hst4.BinEdges(2:end),len-hst4.Values,'color',colors(4,:),'LineWidth',1.5,...
    'DisplayName','3D');
% h(end+1)=plot(hst41.BinEdges(2:end),len-hst41.Values,':','color',colors(4,:),'LineWidth',1.5,...
%     'DisplayName','3D');

%plot(log10([5 5]),[0 len],'k--','LineWidth',1.5);
%plot(log10(sqrt(2/3).*[3 3]),[0 len],'k:','LineWidth',1.5);
xlim([0 20])
xlabelmine('nearest neigbor distance / max($R_\mathrm{vir}$)');
ylabelmine('number of objects with value greater than');
legend(h,'Interpreter','latex')
grid
myAxis;
%plot([3 3],[0 len],'k--','LineWidth',1.5);
%% compare isolatation criteria 
len=length(nneib3.distance);
myFigure;
i=1;
histogram(log10(nneib2.distance(i,:)./rnorm),'Normalization','cumcount');
hold on
i=2;
histogram(log10(nneib2.distance(i,:)./rnorm),'Normalization','cumcount');
i=3;
histogram(log10(nneib2.distance(i,:)./rnorm),'Normalization','cumcount');
histogram(log10(nneib3.distance./max(R200c(dwarfIndx),R200c(allIndx(nneib3.indx)))),'Normalization','cumcount');
%%


