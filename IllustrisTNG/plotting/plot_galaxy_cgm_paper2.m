%% plot images for galaxies
simName='100';
snap=99;
bp=illustris.set_env(simName);

global simDisplayName
global illUnits
global DEFAULT_PRINTOUT_DIR

printoutDirBase=[DEFAULT_PRINTOUT_DIR '/maps'];
if perlimFlag
    loadFofSubTNG100
    
    %      fofs=illustris.groupcat.loadHalos(bp,snap);
    %      subs=illustris.groupcat.loadSubhalos(bp,snap);
    subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
    %
    %     if readFlag
    %         fprintf(' *** Reading data *** \n');
    %         fofs=illustris.groupcat.loadHalos(bp,snap);
    %         subs=illustris.groupcat.loadSubhalos(bp,snap);
    %     end
    %
    %     massThresh=10^9; % threshold for *stellar* mass
    %     ssfr=illustris.utils.calc_ssfr(subs);
    %
    %     massAllGals= double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit); % stellar mass within 2*rhalf
    %     galMask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'centrals','hasGas');
end


idGal=595189 ;
proj='xz';

idHost=subsInfo.hostFof(idGal+1);

printoutDir=[printoutDirBase '/id' num2str(idGal)];
printTag=['id' num2str(idGal) '_Pap'];

%baseName=['galMap_id%i_' proj '%s_snp' num2str(snap) '-' simDisplayName];

galMass=double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,idGal+1).*illUnits.massUnit); % stellar mass within 2*rhalf
hostMass=fofs.Group_M_Crit200(idHost+1).*illUnits.massUnit;

center=subs.SubhaloPos(:,idGal+1);

rhalfGas=subs.SubhaloHalfmassRadType(illustris.partTypeNum('gas')+1,idGal+1);
rhalfStars=subs.SubhaloHalfmassRadType(illustris.partTypeNum('star')+1,idGal+1);
r200=fofs.Group_R_Crit200(idHost+1);
r500=fofs.Group_R_Crit500(idHost+1);

rLimitOuter=1.5*rhalfGas;
rLim=1.05.*rLimitOuter;
boxx=2*rLim;
scales=[5 10 20 30 50 100 200 500];
scaleLength=scales(find(boxx/6>scales,1,'last'));



circ(1).radius=2.*rhalfStars;
circ(1).width=1.5;
circ(1).color=[0 0 0];
circ(1).type='--';
circ(2).radius=rhalfGas;
circ(2).width=1.5;
circ(2).color=[0 0 0];
circ(2).type='-';
% circ(3).radius=r200;
% circ(3).width=1.5;
% circ(3).color=[1 1 1];
% circ(3).type=':';
% circ(4).radius=r200;
% circ(4).width=1.5;



%% load  and  plot star data
if readFlag
    % load stars
    starFields={'Coordinates','Masses','GFM_StellarFormationTime'};
    stars=illustris.snapshot.loadHalo(bp,snap,idHost,'stars',starFields);
    stars.newCoord=illustris.utils.centerObject(stars.Coordinates,center);
    starMask=true(size(stars.newCoord(1,:)));    %   abs(stars.newCoord(1,:))<=rLim & abs(stars.newCoord(2,:))<=rLim & abs(stars.newCoord(3,:))<=rLim ;
    starMask=starMask & stars.GFM_StellarFormationTime>=0; % remove wind particles
    
    
    % load gas data
    gasFields={'Coordinates','Masses','Density','ElectronAbundance'...
        'EnergyDissipation','GFM_CoolingRate','GFM_Metallicity',...
        'InternalEnergy','Machnumber','MagneticField','NeutralHydrogenAbundance',...
        'Potential','StarFormationRate','Velocities'};
    
    gas=illustris.snapshot.loadHalo(bp,snap,idHost,'gas',gasFields);
    
    gas.newCoord=illustris.utils.centerObject(gas.Coordinates,center);
    
    gas=illustris.utils.addEntropy(gas);
    gas=illustris.utils.addTemperature(gas);
    gas=illustris.utils.addPressure(gas);
    
    gasMask=true(size(gas.newCoord(1,:)));
    %gasMask=abs(gas.newCoord(1,:))<=rLim & abs(gas.newCoord(2,:))<=rLim & abs(gas.newCoord(3,:))<=rLim ;
    
end
%vcm=fofs.GroupVel(:,idHost+1); %       sum(gas.Velocities(:,gasMask).*gas.Masses(gasMask),2)./sum(gas.Masses(gasMask));

%% plot
ng=512;
%titleTag=['gal ID: ' num2str(idGal)];

mgalText=mk_exponent_string(galMass);
mhostText=mk_exponent_string(hostMass);

%titleTag=['$' mgalText ',\,'  mhostText '\,\mathrm{M_\odot}$'];
msTag=['$M_\mathrm{gal}=' mgalText '\,\mathrm{M_\odot}$'];
mhTag=['$M_\mathrm{host}=' mhostText '\,\mathrm{M_\odot}$'];


cmapTemp=brewermap(256,'*Spectral');


%cmapTc=brewermap(256,'YlOrRd');
%cmapPress=brewermap(256,'YlOrBr');

caxStar=[4.5 8.5];
caxTemp=[4 6.2];
caxEnt=[-1.2 1.8];
caxDens=[-4.5 -1];
caxTc=[-1.3 1.16];
%caxPress=[-5 -2];
caxTcTff=1.5*[-1 1];
thickness=2*rhalfGas;

scaleArw.start=[-0.8*boxx/2 -0.8*boxx/2];
scaleArw.stop=[-0.8*boxx/2+scaleLength -0.8*boxx/2];
scaleArw.width=2;
scaleArw.color='k';
scaleArw.ends='none';

textS(1).str=['$' num2str(scaleLength) '\,\mathrm{kpc}$'];
textS(1).position=[-0.8*boxx/2 -0.7*boxx/2 ];
textS(1).color=[0 0 0];
textS(1).Rotation=0;

textS(2).str='$z=0$';
%textS(2).position=[0.5*boxx/2 -0.85*boxx/2 ];
textS(2).position=[-0.9*boxx/2 0.82*boxx/2 ];
textS(2).color=[0 0 0];
textS(2).Rotation=0;

textS(3).str=msTag;%   ['$' mgalText '\,\mathrm{M_\odot}$'];
textS(3).position=[-0.9*boxx/2 0.82*boxx/2 ];
textS(3).color=[0 0 0];
textS(3).fontsize=18;
textS(3).Rotation=0;

textS(4).str=['ID: ' num2str(idGal) ', ' simDisplayName];
textS(4).position=[-0.9*boxx/2 0.82*boxx/2 ];
textS(4).color=[0.3 0.3 0.3];
textS(4).fontsize=14;
textS(4).Rotation=0;

textS(6).str=mhTag;%   ['$' mgalText '\,\mathrm{M_\odot}$'];
textS(6).position=[-0.9*boxx/2 0.82*boxx/2 ];
textS(6).color='w';
textS(6).fontsize=18;
textS(6).Rotation=0;


textS(7)=textS(4);
textS(7).Rotation=270;
textS(7).position=[-0.9*boxx/2 0.5*boxx/2 ];
textS(7).fontsize=14;

barProp.fontsize=14;
%barProp.location='east';
barProp.position=[0.785,0.135,0.028,0.765];
barProp.color='k';
barProp.tickLabelInterpreter='latex';
barProp.axisLocation='in';
%% stars

textS(5).str='$\log \rho_\mathrm{star}\,[\mathrm{M_\odot/kpc^2}]$';
textS(5).position=[-0.43*boxx/2 -0.87*boxx/2 ];
textS(5).color='k';
textS(5).fontsize=18;
textS(5).Rotation=0;


cmapStar=brewermap(256,'Greys');

illustris.plots.mkmapStars('star',stars,'type','mass','ng',512,proj,'nn',3,...
    'mask',starMask,'box',boxx,'cmap',cmapStar,'white','labels','no','arrow',scaleArw,'text',textS([1 3 5 ]),...
    'zoom',rLimitOuter,'clims',caxStar,'circ',circ,'noticks','nogrid','barprop',barProp,'bartag',' ','print',printTag,'savefig','printoutdir',printoutDir);

%% gas


%% temp
%caxTemp=[4 6.5];
textS(5).position=[0.15*boxx/2 -0.87*boxx/2 ];
textS(5).str='$\log T\,[\mathrm{K}]$';

illustris.plots.mkmapGas('gas',gas,'type','temp','ng',ng,proj,...
    'mask',gasMask,'circ',circ,'clims',caxTemp,'cmap',cmapTemp,'nanval','max',...%'vfield','dilute',6,'vcm',vcm,...
    'box',boxx,'zoom',rLimitOuter,'thick',thickness,'text',textS([1 4 5]),...
    'arrow',scaleArw,'noticks','labels','no','nogrid','barprop',barProp,'bartag',' ','print',printTag,'savefig','printoutdir',printoutDir);

%ht=titlemine('Temperature',16);set(ht,'color','w');

%% tcool 2
textS(5).position=[0.01*boxx/2 -0.87*boxx/2 ];
textS(5).str='$\log t_\mathrm{cool}/t_\mathrm{ff}$';
textS(3).color='k';
%textS(5).color='w';
%barProp.color='w';

cmapTcTff=brewermap(256,'*RdBu');
illustris.plots.mkmapGas('gas',gas,'id',idGal,'r200',r200,'type','tctff','ng',ng,proj,...
    'mask',gasMask,'circ',circ,'clims',caxTcTff,'cmap',cmapTcTff,'nanval','max',...%'vfield','dilute',6,'vcm',vcm,...
    'box',boxx,'zoom',rLimitOuter,'thick',thickness,'nobackground',...
    'text',textS([1 3 5 7]),'arrow',scaleArw,'noticks','labels','no','nogrid','barprop',barProp,'bartag',' ','print',printTag,'savefig','printoutdir',printoutDir);
%'title',titleTag,'print',printTag,'savefig','printoutdir',printoutDir); %,...%'streamdense',0,...

%% pressure

% cmapPress=brewermap(256,'YlOrRd');
% %brewermap(256,'YlOrBr');
%
% illustris.plots.mkmapGas('gas',gas,'type','pressure','ng',ng,proj,...
%     'mask',gasMask,'circ',circ,'cmap',cmapPress,'clims',caxPress,...%'vfield','dilute',6,'vcm',vcm,...
%     'box',boxx,'zoom',rLimitOuter,'thick',thickness,'nobackground','text',textS,'arrow',scaleArw,'noticks','labels','no','nogrid',...
%     'print',printTag,'savefig','printoutdir',printoutDir);
%

textS(1).color=[1 1 1];
textS(2).color=[1 1 1];
textS(3).color=[1 1 1];
scaleArw.color='w';
barProp.color='w';
%% density
cmapDens=brewermap(256,'*YlGnBu');
%cmapDens=viridis(256);
textS(5).position=[-0.05*boxx/2 -0.87*boxx/2 ];
textS(5).str='$\log n\,[\mathrm{cm^{-3}}]$';
textS(5).color='w';

illustris.plots.mkmapGas('gas',gas,'type','ndensity','ng',ng,proj,...
    'mask',gasMask,'circ',circ,'clims',caxDens,'cmap',cmapDens,...%'vfield','dilute',6,'vcm',vcm,...
    'box',boxx,'zoom',rLimitOuter,'thick',thickness,'text',textS([1 2 5]),'arrow',scaleArw,...
    'noticks','labels','no','nogrid','barprop',barProp,'bartag',' ','print',printTag,'savefig','printoutdir',printoutDir);


%% entropy
%cmapEnt=brewermap(256,'PuRd');
cmapEnt=flipud(plasma(256));
textS(5).position=[-0.1*boxx/2 -0.87*boxx/2 ];
textS(5).str='$ S\,[\mathrm{KeV\, cm^2}]$';

illustris.plots.mkmapGas('gas',gas,'type','entropy','ng',ng,proj,...
    'mask',gasMask,'circ',circ,'clims',caxEnt,'cmap',cmapEnt,'nanval','max',...%'vcm',vcm,'streamdense',1.5,...    %'vfield','dilute',10,'vcm',vcm,...
    'box',boxx,'zoom',rLimitOuter,'thick',thickness,'text',textS([1 5 6]),'arrow',scaleArw,...
    'noticks','labels','no','nogrid','barprop',barProp,'bartag',' ','print',printTag,'savefig','printoutdir',printoutDir);

%% tcool
cmapTc=flipud(viridis(256)); %brewermap(256,'*YlOrRd');
%cmapTc=brewermap(256,'YlGnBu');  %brewermap(256,'YlOrRd');
textS(5).position=[-0.1*boxx/2 -0.87*boxx/2 ];
textS(5).str='$\log t_\mathrm{cool}\,[\mathrm{Gyr}]$' ;

illustris.plots.mkmapGas('gas',gas,'type','tcool','ng',ng,proj,...
    'mask',gasMask,'circ',circ,'clims',caxTc,'cmap',cmapTc,...
    'box',boxx,'zoom',rLimitOuter,'thick',thickness,'nobackground','nanval','max',...
    'text',textS([1 5]) ,'arrow',scaleArw,'noticks','labels','no','nogrid',...
    'barprop',barProp,'bartag',' ','print',printTag,'savefig','printoutdir',printoutDir);








