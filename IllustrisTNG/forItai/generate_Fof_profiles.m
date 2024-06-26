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
m200=fofs.Group_M_Crit200(fofID+1);
r200=fofs.Group_R_Crit200(fofID+1);
center=fofs.GroupPos(:,fofID+1);
rFactor=2.5;


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
%tic;
sumPart=0;
for i=1:nChunks
    %fprintf('i=%i \n',i);
    % build the subset structure     
    startPoint=(i-1).*chunkLength;
    endPoint=min(i.*chunkLength-1,numPartTotal);
    subset.offsetType(nPartType)=startPoint;
    subset.lenType(nPartType)=endPoint-startPoint;

    % read in the gas particles in the chunk 
    gas=illustris.snapshot.loadSubset(BASEPATH, snap,'gas',gasFields,subset);

    %fprintf('loaded %i gas cells \n',gas.count);
    % idenitfy the relevant particles
    
    gas.newCoord = illustris.utils.centerObject(gas.Coordinates,center);
    gasDist=sqrt( sum(double(gas.newCoord).^2,1));
    mask=gasDist<=rmax;
    sumPart=sumPart+sum(mask);
    if any(mask)
        gas=mask_structure(gas,mask);
        gas.count=length(gas.Masses);
        fprintf('Bingo! %i, %i %i \n',i,sum(mask),gas.count)
        % add to list
        if ~firstFlag
            gasCells=illustris.infrastructure.concat_particle_struct(gasCells,gas);
        else
            gasCells=gas;
            firstFlag=false;
        end
    end
    %toc;
    %fprintf('press\n');pause;
end

fprintf('sumPart = %i, count= %i \n',sumPart,gasCells.count);


%% calculate profiles

% set radial profile 
nbins=1000;
rp=linspace(0,rmax,nbins);

% add temperature
gasCells=illustris.utils.addTemperature(gasCells);

% set cell radius 
gasCells=illustris.utils.addCellRadius(gasCells);

%gasDist=sqrt( sum(double(gasCells.newCoord).^2,1));
vol=gasCells.Masses./gasCells.Density.*illUnits.volumeUnit;
mass=gasCells.Masses.*illUnits.massUnit;
coords=gasCells.newCoord.*illUnits.lengthUnit;
cellRad=gasCells.CellRadius.*illUnits.lengthUnit;
profs.density = mk_radial_profile_cells(coords,cellRad,...
    gasCells.Density.*illUnits.densityUnit,...
    'bins',rp,'wt',vol,'type','intensiveFancy');
profs.temperature = mk_radial_profile_cells(coords,cellRad,...
    gasCells.Temperature,...
    'bins',rp,'wt',mass,'type','intensiveFancy');




%% save to mat files. 
fprintf('save?\n');pause;
global DRACOFLAG
if DRACOFLAG
    global DEFAULT_MATFILE_DIR
    global simDisplayName
    %fname=sprintf('gasProperties_z%s_%s',num2str(illustris.utils.get_zred(snap)),simDisplayName);
    fname=sprintf('gasCutoutTest_snp%s_%s.mat',num2str(snap),simDisplayName);
    save([DEFAULT_MATFILE_DIR '/' fname],'gasCells','profs','-v7.3')
    
    fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
    
    
end

%% generate profiles  

% set radial profile 

%gasCells=illustris.utils.addTemperature(gasCells);

% set cell radius 

%gas.newCoord = illustris.utils.centerObject(gas.Coordinates,center);
%gasDist=sqrt( sum(double(gas.newCoord).^2,1));

 %res = mk_radial_profile_cells(cellPos,cellRad,val,varargin)
