bp=illustris.set_env(50);

snap=67
idGal=477099

loadFofSub

global illUnits
illustris.utils.set_illUnits(snap);
illUnits

center=subs.SubhaloPos(:,idGal+1);
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
idHost=subsInfo.hostFof(idGal+1)
rhalfGas=subs.SubhaloHalfmassRadType(illustris.partTypeNum('gas')+1,idGal+1);
rhalfStars=subs.SubhaloHalfmassRadType(illustris.partTypeNum('star')+1,idGal+1);
r200=fofs.Group_R_Crit200(idHost+1);
r500=fofs.Group_R_Crit500(idHost+1);
m200=fofs.Group_M_Crit200(idHost+1).*illUnits.massUnit
% gas=illustris.snapshot.loadHalo(bp,snap,idHost,'gas',gasFields);
% gas.newCoord=illustris.utils.centerObject(gas.Coordinates,center);
gasFields={'Coordinates','Masses','Density','ElectronAbundance'...
'EnergyDissipation','GFM_CoolingRate','GFM_Metallicity',...
'InternalEnergy','Machnumber','MagneticField','NeutralHydrogenAbundance',...
'Potential','StarFormationRate','Velocities'};
gas=illustris.snapshot.loadHalo(bp,snap,idHost,'gas',gasFields);
%gas=illustris.snapshot.loadSubhalo(bp,snap,idGal+1,'gas',gasFields);
gas.newCoord=illustris.utils.centerObject(gas.Coordinates,center);
gas=illustris.utils.addEntropy(gas);
gas=illustris.utils.addTemperature(gas);
gas=illustris.utils.addPressure(gas);
ng=512;

%rhalfGas.*illUnits.lengthUnit
boxx=4*rhalfGas;

illustris.plots.mkmapGas('gas',gas,'type','ndensity','ng',ng,'thick',-1*ng,...
    'print',['id' num2str(idGal)],'cmap',inferno(256))