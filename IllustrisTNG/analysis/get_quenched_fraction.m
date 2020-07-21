
% load stuff
sim=100;
bp=illustris.set_env(sim);

global illUnits 

load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/fofs_subs_TNG100_z0.mat')

subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs,snap);

massThresh=1e9;
mask=illustris.utils.generateMask('subs',subs,'fofs',fofs,'sats','mass',massThresh);

% get ssfr 
ssfr=illustris.utils.calc_ssfr(subs);
ssfr=ssfr(mask);

% get mass 
galMass=subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,mask).*illUnits.massUnit;

% get host mass 
hostMass=fofs.Group_M_Crit200(subs.SubhaloGrNr(mask)+1).*illUnits.massUnit;

% get rad
pos=subs.SubhaloPos(:,mask)-fofs.GroupPos(:,subs.SubhaloGrNr(mask)+1);
galRad=sqrt(sum(pos.^2,1))./fofs.Group_R_Crit200(subs.SubhaloGrNr(mask)+1);


% define quenched 
quenchedThresh=1e-11;
qMask=ssfr<quenchedThresh;



