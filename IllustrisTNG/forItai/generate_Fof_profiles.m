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

%% identify the gas cells we are interested in for the halo 
% load chunks of the snapshot data, find the relevant particles and create
% a subset list to work with later 

gasFields={'Coordinates','Density','ElectronAbundance','InternalEnergy','Masses'};%,'Velocities','StarFormationRate'};

% load the header

% identify the no. of particles 



% build the subset structure 
subset=...
% load a chunk of particles 

gas=illustris.snapshot.loadSubset(BASEPATH, snap,'gas',gasFields,subset);
% idenitfy the relevant particles 

gas.newCoord = illustris.utils.centerObject(gas.Coordinates,center);
gasDist=sqrt( sum(double(gas.newCoord).^2,1));
mask=gasDist<=rmax;

gas=mask_structure(ga,mask);

% add to list 
gasCells=illustris.infrastructure.concat_particle_struct(gasCells,gas2);






%% generate profiles  

% set radial profile 

gasCells=illustris.utils.addTemperature(gasCells);

% set cell radius 

%gas.newCoord = illustris.utils.centerObject(gas.Coordinates,center);
%gasDist=sqrt( sum(double(gas.newCoord).^2,1));

 res = mk_radial_profile_cells(cellPos,cellRad,val,varargin)
