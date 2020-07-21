
bp=illustris.set_env(100);


%% load  BH stuff. 
snap=99;


illustris.utils.set_illUnits(snap);

global illUnits
global DEFAULT_MATFILE_DIR
global simDisplayName

loadFofSubTNG100

subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
units; %

load([DEFAULT_MATFILE_DIR '/BH_energyInjection_snp' num2str(snap) '_' simDisplayName '.mat'])
    


%% load ID's 

loadClusterIDs;
loadControlIDs;

clusterIndx=clustergalaxyids+1;
controlIndx=controlgalaxyids+1;

cluster.stellarMass=bhStruct.galMass(clusterIndx);
cluster.bhQM=bhStruct.inGal.cumEngQM(clusterIndx);
cluster.bhRM=bhStruct.inGal.cumEngRM(clusterIndx);
cluster.bhMass=bhStruct.inGal.bhMassMax(clusterIndx);

control.stellarMass=bhStruct.galMass(controlIndx);
control.bhQM=bhStruct.inGal.cumEngQM(controlIndx);
control.bhRM=bhStruct.inGal.cumEngRM(controlIndx);
control.bhMass=bhStruct.inGal.bhMassMax(controlIndx);


%% write to file 
outDir='/home/zinger/workProjects/matlab_scripts/IllustrisTNG/forGandhali/';
fnameCluster=[outDir 'BHproperties_cluster_snp' num2str(snap) '_' simDisplayName '.txt'];
fnameControl=[outDir 'BHproperties_control_snp' num2str(snap) '_' simDisplayName '.txt'];

fidClust=fopen(fnameCluster,'w');
fidCntrl=fopen(fnameControl,'w');
header{1}='## stellar mass, BH mass, and energy injection all within 2 times stellar half mass radius';
header{3}='## Columns: Gal ID, Stellar mass, BH mass, E_QM, E_RM';
header{2}='## Units: mass = solar mass, energy: 1e53 ergs';

% header
for i=1:3
    fprintf(fidClust,'%s \n',header{i});
    fprintf(fidCntrl,'%s \n',header{i});
end

% data 
format='%i   %6.4g   %6.4g   %6.4g   %6.4g \n';
for i=1:length(clustergalaxyids)
    fprintf(fidClust,format,clustergalaxyids(i),cluster.stellarMass(i),cluster.bhMass(i),...
        cluster.bhQM(i),cluster.bhRM(i));
end

for i=1:length(controlgalaxyids)
    fprintf(fidCntrl,format,controlgalaxyids(i),control.stellarMass(i),control.bhMass(i),...
        control.bhQM(i),control.bhRM(i));
end

fclose(fidClust);
fclose(fidCntrl);


    



