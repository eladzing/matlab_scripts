%% generate a catalog for hydrogen gas in the different regions of subahlos. 


illustris.utils.set_illUnits(snap);

global illUnits
global DEFAULT_MATFILE_DIR
global simDisplayName

if ~exist('readFlag','var')
    readFlag=true;
end

%% read catalogs
if readFlag
    fprintf(' *** Reading data *** \n');
     
    % NEED TO ADD FIELDS!
    subFields={'SubhaloMassInRadType','SubhaloFlag','SubhaloGrNr',...
        'SubhaloMassType','SubhaloMassInRad', 
    
    fofFields={'GroupFirstSub','Group_M_Crit500','Group_M_Crit200',...
    'Group_M_Mean200','Group_M_TopHat200','Group_R_Crit200'
    
    
    fofs=illustris.groupcat.loadHalos(bp,snap);
    subs=illustris.groupcat.loadSubhalos(bp,snap);
    
    readFlag=false;
    
    gasFields={'Coordinates','ParticleIDs'
    
end

subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
units; % load general unit structure in cgs.

%% select galaxies 
massAllGals= double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit); % stellar mass within 2*rhalf

% this mask selects all galaxies with dm component & stars, above *stellar* mass limit whose host has virials parameters
galMask=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap);


%% select radii 
r200c=fofs.Group_R_Crit200(subsInfo.hostFof+1);

%% read hydrogen catalog 

hydrogenMass=illustris.utils.read_hydrogen_catalog(snap);
%% run over halos
cnt=0;
for id=0:len-1
 perCent=floor((id+1)/len*100);
    
    if (perCent>=stepNext)     %       mod(perCent,5)==0 && perCent>=5)
        fprintf('%s %% of Galaxies done  \n',num2str(floor(perCent)))
        stepNext=stepNext+step;
    end
    
    if galMask(id+1)
        
        
        cnt=cnt+1;
        gas=illustris.snapshot.loadSubhalo(bp, snap, id, 'gas',gasFields);
%% sum up halos 