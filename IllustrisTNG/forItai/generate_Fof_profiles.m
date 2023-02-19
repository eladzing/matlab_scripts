%% This script generates 1D profiles of density and temperature 
% out to several times R200 of a given Fof.

%% halo ID and parameters 
%sim='TNG50';
snap=99;
fofID=5;

%% load relevant environmet 
bp=illustris.set_env(103);
illustris.utils.set_illUnits(snap) % set the simulation units for the snapshot.

[subs,fofs,subsInfo]=illustris.loadFofSub(snap);

global illUnits
global cosmoStruct
global BASEPATH

%% get relavent Fof info 
r200=fofs.Group_M_Crit200(fofID+1);
m200=fofs.Group_R_Crit200(fofID+1);
center=fofs.GroupPos(:,fofID+1);
rFactor=1;


rmax=r200.*rFactor;

%% identify the gas cells we are interested in for the halo 
% load chunks of the snapshot data, find the relevant particles and create
% a subset list to work with later 
nPartType=illustris.partTypeNum('gas')+1;
gasFields={'Coordinates','Density','ElectronAbundance','InternalEnergy','Masses'};%,'Velocities','StarFormationRate'};

nChunks=20; % no. of chumks to break up the particle acquisition


% load the header and get the number of particles 
 header = illustris.snapshot.loadHeader(BASEPATH,snap);
 nPart = illustris.snapshot.getNumPart(header); % no. of particles in the entire sim 
numPartTotal=nPart(nPartType);
chunkLength=numPartTotal./nChunks;

%get the snap offsets for the subset
subset= illustris.snapshot.getSnapOffsets(BASEPATH,snap,0,'Subhalo'); 



%% load chunks of particles 
firstFlag=true;
tic;
for i=1:nChunks
    fprintf('i=%i \n',i);
    % build the subset structure     
    startPoint=(i-1).*chunkLength
    endPoint=min(i.*chunkLength-1,numPartTotal)
    subset.offsetType(nPartType)=startPoint;
    subset.lenType(nPartType)=endPoint-startPoint

    % read in the gas particles in the chunk 
    gas=illustris.snapshot.loadSubset(BASEPATH, snap,'gas',gasFields,subset);

    fprintf('loaded %i gas cells \n',gas.count);
    % idenitfy the relevant particles
    
    gas.newCoord = illustris.utils.centerObject(gas.Coordinates,center);
    gasDist=sqrt( sum(double(gas.newCoord).^2,1));
    mask=gasDist<=rmax;
    
    if sum(mask)>0
        gas=mask_structure(gas,mask);
        fprintf('Bingo!')
        % add to list
        if ~firstFlag
            gasCells=illustris.infrastructure.concat_particle_struct(gasCells,gas);
        else
            gasCells=gas;
            firstFlag=false;
        end
    end
    toc;
end

%% save to mat files. 

global DRACOFLAG
if DRACOFLAG
    global DEFAULT_MATFILE_DIR
    global simDisplayName
    %fname=sprintf('gasProperties_z%s_%s',num2str(illustris.utils.get_zred(snap)),simDisplayName);
    fname=sprintf('gasCutoutTest_snp%s_%s.mat',num2str(snap),simDisplayName);
    save([DEFAULT_MATFILE_DIR '/' fname],'tCoolStruct','-v7.3')
    
    fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
    
    
end

%% generate profiles  

% set radial profile 

%gasCells=illustris.utils.addTemperature(gasCells);

% set cell radius 

%gas.newCoord = illustris.utils.centerObject(gas.Coordinates,center);
%gasDist=sqrt( sum(double(gas.newCoord).^2,1));

 %res = mk_radial_profile_cells(cellPos,cellRad,val,varargin)
