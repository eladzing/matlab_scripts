%% This script generates 1D profiles of density and temperature 
% out to several times R200 of a given Fof.

%% halo ID and parameters 
sim='TNG50';
snap='';
fofID=0;

%% load relevant environmet 
bp=illustris.set_env(100);
illustris.utils.set_illUnits(snap) % set the simulation units for the snapshot.

[subs,fofs,subsInfo]=illustris.loadFofSub(snap);

global illUnits
global cosmoStruct
global BASEPATH

%% get relavent Fof info 
r200=fofs.Group_M_Crit200(fofID);
m200=fofs.Group_R_Crit200(fofID);
center=fofs.GroupPos(:,fofID);
rFactor=5;


rmax=r200.*rFactor;

%% load gas cells 
gasFields={'Coordinates','Density','ElectronAbundance','InternalEnergy','Masses'};%,'Velocities','StarFormationRate'};

gas=illustris.snapshot.loadSubset(BASEPATH, snap,'gas',gasFields);

gas=illustris.utils.addTemperature(gas);

gas.newCoord = illustris.utils.centerObject(gas.Coordinates,center);
gasDist=sqrt( sum(double(gas.newCoord).^2,1));

%% identify all cells within the maximal radius 
mask=gasDist




