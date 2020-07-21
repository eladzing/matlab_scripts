%% plot single galaxy 
bp=illustris.set_env('100');
subID=143887;
fofID=22;
snap=99;
load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/tng100_z0_fofs_subs.mat')

gas=illustris.snapshot.loadSubhalo(bp,snap,subID, 'gas');
stars=illustris.snapshot.loadSubhalo(bp,snap,subID, 'stars');
%subs=illustris.groupcat.loadSubhalos(bp,snap);fofs=illustris.groupcat.loadHalos(bp,snap);

gas.newCoord=illustris.utils.centerObject(gas.Coordinates,subs.SubhaloPos(:,subID+1));
stars.newCoord=illustris.utils.centerObject(stars.Coordinates,subs.SubhaloPos(:,subID+1));
bx=10.*subs.SubhaloHalfmassRadType(1,subID+1);

illustris.plots.mkmapGas('gas',gas,'ng',512,'box',