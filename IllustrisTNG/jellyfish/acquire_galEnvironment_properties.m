%% get the environmental conditions around the galaxies 

global illUnits
global DEFAULT_MATFILE_DIR

load([DEFAULT_MATFILE_DIR '/cosmic_jellyfish_objectTable.mat'])
snaps=unique(objectTable.snap);
sims=unique(objectTable.sim);

%% set aperatures 

aper=[];


%%  set FOF gas fields, 
fields={'Masses','Density','Coordinates','NeutralHydrogenAbundance','StarFormationRate',...
    'ElectronAbundance','EnergyDissipation','InternalEnergy',...
    'Machnumber','Velocities'}; %% add metalicity


%% initialize

% snap & sim and id
galFoFProps.snap=objectTable.snap;
galFoFProps.sim=objectTable.sim;
galFoFProps.subfind=objectTable.subfind;

ngal=length(objectTable.subfind);

galFoFProps.hostID=zeros(1,ngal);
galFoFProps.nSatHost=zeros(1,ngal); % no of satellites in host
galFoFProps.rhoLocal=zeros(length(aper),size(objectTable.subfind));
galFoFProps.ramPress=zeros(length(aper),size(objectTable.subfind));
galFoFProps.mach=zeros(length(aper),size(objectTable.subfind));
galFoFProps.velMedium=zeros(length(aper),3,size(objectTable.subfind));
galFoFProps.csonic=zeros(length(aper),size(objectTable.subfind));

%%run over simulations and snapshots
for k=1:length(sims)
    
    bp=illustris.set_env(sims(k));
    %global simDisplayName
    global LBox
    fprintf('running on simulation %s \n ', sims(k));
    
    simMask=strcmp(objectTable.sim,sims(k));
    
    for i=1:length(snaps)
        
        snap=snaps(i);
        
        fprintf('Loading catalogs of snap %i \n',snap);
        
        [subs,fofs,subsInfo]=illustris.loadFofSub(snap);
        
        illustris.utils.set_illUnits(snap);
                
        fprintf('Getting galaxy properties \n');
        
        % find indices of gals from the snapshot
        snapMask=objectTable.snap==snap;
        %fullMask=snapMask & simMask;
        inds=find(snapMask & simMask);
        
        galInds=objectTable.subfind(inds)+1; % indices in the catalogs 
        hostInds=subsInfo.hostFof(galInds)+1;% indices in the catalogs 
        
        hostList=unique(hostInds); % list of hosts in this particular snapshot/sim
        
        %% loop over hosts 
        for j=1:length(hostList)
            
            hind=hostList(j);
            % find the relavant sat's in this particular host 
            
            
            % load FOF gas for a host 
            gas=illustris.snapshot.loadHalo(bp,snap,hind,'gas',fields);
        
            % calculate temperature & entropy 
            gas=illustris.utils.addTemperature(gas);
            gas=illustris.utils.addEntropy( gas );
            
            % center coordinates around host position
            gas.newCoord = illustris.utils.centerObject(gas.CenterOfMass,fofs.GroupPos(:,hind));
            
            
            % fix velocities
            hostVel=fofs.GroupVel(:,hind).*illustris.utils.velocityFactor(snap,'host');
            gas.Velocities=gas.Velocities.*illustris.utils.velocityFactor(snap,'gas');
            for kk=1:3
                gas.Velocities(kk,:)=gas.Velocities(kk,:)-hostVel(kk);
            end
            
            
            
            
        %% fill out the properties
        
        
        
        % virial properties
        galFoFProps.hostID(inds)=subsInfo.hostFof(galInds);
        
        galFoFProps.hostM200c(inds)=illUnits.massUnit.*...
             double(fofs.Group_M_Crit200(hostInds));
        
        galFoFProps.hostR200c(inds)=illUnits.lengthUnit.*...
             double(fofs.Group_R_Crit200(hostInds));
        
        galFoFProps.hostM200m(inds)=illUnits.massUnit.*...
             double(fofs.Group_M_Mean200(hostInds));
        
        galFoFProps.hostR200m(inds)=illUnits.lengthUnit.*...
             double(fofs.Group_R_Mean200(hostInds));
        
        % gal/subfind properties
        galFoFProps.stellarMass(inds)=illUnits.massUnit.*...
             double(subs.SubhaloMassType(illustris.partTypeNum('stars')+1,galInds));
        
        galFoFProps.galStellarMass(inds)=illUnits.massUnit.*...
             double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,galInds));
        
        galFoFProps.gasMass(inds)=illUnits.massUnit.*...
             double(subs.SubhaloMassType(illustris.partTypeNum('gas')+1,galInds));
        
        galFoFProps.galGasMass(inds)=illUnits.massUnit.*...
             double(subs.SubhaloMassInRadType(illustris.partTypeNum('gas')+1,galInds));
        
        galFoFProps.galBHMass(inds)=illUnits.massUnit.*...
             double(subs.SubhaloMassInRadType(illustris.partTypeNum('BH')+1,galInds));
        
        galFoFProps.galSFR(inds)=...
             double(subs.SubhaloSFRinRad(galInds));
        
       %% position data 
       galPos=double(subs.SubhaloPos(:,galInds)); % global position in simulation box, in simulation units 
       hostPos=double(fofs.GroupPos(:,hostInds));% global position of host in simulation box, in simulation units 
       
       
       % apply periodic boundary conditions to calculate the actual radial
       % position 
       galFoFProps.pos(:,inds)=galPos; 
       
       rpos=findDistance(galPos,hostPos,LBox,3);
       galFoFProps.rpos(inds)=rpos.*illUnits.lengthUnit;
       
       
       
       % leave this for now till we test.
       for ii=1:3
          mask=abs(galPos(ii,:)-hostPos(ii,:))>0.5.*LBox;
          
          if any(mask)
            mask1=mask & hostPos(ii,:)>0.5*LBox;
            mask2=mask & hostPos(ii,:)<=0.5*LBox;
            
            galPos(ii,mask1)=galPos(ii,mask1)+LBox;
            galPos(ii,mask2)=galPos(ii,mask2)-LBox;
          end
          
          galPos(ii,:)=galPos(ii,:)-hostPos(ii,:);
             
       end
       
       % velocity w.r.t host
       vsat=double(subs.SubhaloVel(:,galInds)).*illustris.utils.velocityFactor(snap,'sub'); % 3d velocity for each galaxy 
       vhost=double(fofs.GroupVel(:,hostInds)).*illustris.utils.velocityFactor(snap,'host');% 3d velocity of host for each galaxy 
        
       vel=vsat-vhost;
              
       
       
%        vrad=sum(vel.*galPos,1)./rpos;
      
       galFoFProps.vel(:,inds)=vel;
       galFoFProps.vrad(inds)=sum(vel.*galPos,1)./rpos;
       
        
    end
end

fname=sprintf('jf_galProperties_CJF.mat');
save([DEFAULT_MATFILE_DIR '/' fname],'galFoFProps','-v7.3')

fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);



    
    
    % galFoFProps.vel(inds,3)=...
    %
    %     galFoFProps.rpos(inds)=...
    %
    %     galFoFProps.vrad(inds)=...
    %
    %     galFoFProps.vtan(inds)=...
    
    %massThresh=10^8.5; % threshold for *stellar* mass
    
    
    
    
