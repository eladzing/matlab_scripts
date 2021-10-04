%% get the environmental conditions around the galaxies 

global illUnits
global DEFAULT_MATFILE_DIR

load([DEFAULT_MATFILE_DIR '/cosmic_jellyfish_objectTable.mat'])
snaps=unique(objectTable.snap);
sims=unique(objectTable.sim);

%% set aperatures 

aper=[50 100];  % in physical kpc 


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
        tabInds=find(snapMask & simMask);
        
        galInds=objectTable.subfind(tabInds)+1; % indices in the catalogs 
        
        hostInds=subsInfo.hostFof(galInds)+1;% indices in the catalogs 
        
        hostList=unique(hostInds); % list of hosts in this particular snapshot/sim
        
        %% loop over hosts 
        for j=1:length(hostList)
            
            hind=hostList(j);
            % find the relavant sat's in this particular host 
            localGalList=galInds(hostInds==hind);
            localTabInds=tabInds(hostInds==hind);  
            
            
            % load FOF gas for a host 
            gas=illustris.snapshot.loadHalo(bp,snap,hind-1,'gas',fields);
        
            % calculate temperature & entropy 
            gas=illustris.utils.addTemperature(gas);
            gas=illustris.utils.addEntropy( gas );
            
            % center coordinates around host position
            %gas.newCoord = illustris.utils.centerObject(gas.CenterOfMass,fofs.GroupPos(:,hind));
            
            
            % fix velocities and center. 
            hostVel=fofs.GroupVel(:,hind).*illustris.utils.velocityFactor(snap,'host');
            gas.newVel=gas.Velocities.*illustris.utils.velocityFactor(snap,'gas');
            for kk=1:3
                gas.newVel(kk,:)=gas.newVel(kk,:)-hostVel(kk);
            end
            
            %% loop of local galaxies
            for j=1:length(localGalList)
                
                % find index of galaxy in objectTable 
                locInd=localTabInds(j);
                
                %center fof gas around galaxy
                locCoord=illustris.utils.centerObject(gas.CenterOfMass,subs.SubhaloPos(:,localGalList(j)));
                gasDist=sqrt( sum(double(locCoord).^2,1)).*illUnits.lengthUnit;  %distance in physical kpc
                
                galVel=subs.SubhaloVel(:,localGalList(j)).*illustris.utils.velocityFactor(snap,'sub');
                
                %% loop over aperatures
                for ii=1:length(aperatures)
                    
                    aper=apertures(ii);
                    
                    % generate mask around object
                    mask=gasDist<=aper;
                    
                    % local density 
                    rhoLocal=mean(gas.Density(mask));
                    galFoFProps.rhoLocal(ii,locInd)=rhoLocal;
                    
                    % local Velocity
                    for kk=1:3
                        velLocal(kk)=sum(gas.Masses(mask).*gas.newVel(kk,mask));
                    end
                    velLocal=velLocal./sum(gas.Masses(mask));
                    galFoFProps.velMedium(ii,:,locInd)=velLocal;
                    
                    % ram Pressure 
                    galFoFProps.ramPress(ii,locInd)=rhoLocal.*sum((galVel-velLocal).^2);
                    
                    
                    %% fill out the properties
        
            % local density
            
        
        
        
        
        
        
        
        % virial properties
        galFoFProps.hostID(tabInds)=subsInfo.hostFof(galInds);
        
        galFoFProps.hostM200c(tabInds)=illUnits.massUnit.*...
             double(fofs.Group_M_Crit200(hostInds));
        
        galFoFProps.hostR200c(tabInds)=illUnits.lengthUnit.*...
             double(fofs.Group_R_Crit200(hostInds));
        
        galFoFProps.hostM200m(tabInds)=illUnits.massUnit.*...
             double(fofs.Group_M_Mean200(hostInds));
        
        galFoFProps.hostR200m(tabInds)=illUnits.lengthUnit.*...
             double(fofs.Group_R_Mean200(hostInds));
        
        % gal/subfind properties
        galFoFProps.stellarMass(tabInds)=illUnits.massUnit.*...
             double(subs.SubhaloMassType(illustris.partTypeNum('stars')+1,galInds));
        
        galFoFProps.galStellarMass(tabInds)=illUnits.massUnit.*...
             double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,galInds));
        
        galFoFProps.gasMass(tabInds)=illUnits.massUnit.*...
             double(subs.SubhaloMassType(illustris.partTypeNum('gas')+1,galInds));
        
        galFoFProps.galGasMass(tabInds)=illUnits.massUnit.*...
             double(subs.SubhaloMassInRadType(illustris.partTypeNum('gas')+1,galInds));
        
        galFoFProps.galBHMass(tabInds)=illUnits.massUnit.*...
             double(subs.SubhaloMassInRadType(illustris.partTypeNum('BH')+1,galInds));
        
        galFoFProps.galSFR(tabInds)=...
             double(subs.SubhaloSFRinRad(galInds));
        
       %% position data 
       galPos=double(subs.SubhaloPos(:,galInds)); % global position in simulation box, in simulation units 
       hostPos=double(fofs.GroupPos(:,hostInds));% global position of host in simulation box, in simulation units 
       
       
       % apply periodic boundary conditions to calculate the actual radial
       % position 
       galFoFProps.pos(:,tabInds)=galPos; 
       
       rpos=findDistance(galPos,hostPos,LBox,3);
       galFoFProps.rpos(tabInds)=rpos.*illUnits.lengthUnit;
       
       
       
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
      
       galFoFProps.vel(:,tabInds)=vel;
       galFoFProps.vrad(tabInds)=sum(vel.*galPos,1)./rpos;
       
        
    end
end

fname=sprintf('jf_galEnvironmentProperties_CJF.mat');
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
    
    
    
    
