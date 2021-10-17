%% get the environmental conditions around the galaxies 

global illUnits
global DEFAULT_MATFILE_DIR

load([DEFAULT_MATFILE_DIR '/cosmic_jellyfish_objectTable.mat'])
snaps=unique(objectTable.snap);
sims=unique(objectTable.sim);

%% set aperatures 
% positive values in physical kpc, negative values are in units of stellar half mass radius 
aper=[50 100 -20];  


units
%%  set FOF gas fields, 
fields={'Masses','Density','Coordinates','NeutralHydrogenAbundance','StarFormationRate',...
    'ElectronAbundance','EnergyDissipation','InternalEnergy',...
    'Machnumber','Velocities'}; %% add metalicity

% factor for c_sonic calculation
gamma=5/3;
csFac=sqrt(gamma.*Units.kb./Units.muMass./Units.mp)./Units.km;
%% initialize

% snap & sim and id
galFoFProps.snap=objectTable.snap;
galFoFProps.sim=objectTable.sim;
galFoFProps.subfind=objectTable.subfind;

ngal=length(objectTable.subfind);

galFoFProps.hostID=zeros(1,ngal);
galFoFProps.nSatHost=zeros(1,ngal); % no of satellites in host

% full signifies that values are calculated over all relavent gas cells 
galFoFProps.full.rhoLocal=zeros(length(aper),size(objectTable.subfind));
galFoFProps.full.ramPress=zeros(length(aper),size(objectTable.subfind));
galFoFProps.full.mach=zeros(length(aper),size(objectTable.subfind));
galFoFProps.full.velMedium=zeros(length(aper),3,size(objectTable.subfind));
galFoFProps.full.csonic=zeros(length(aper),size(objectTable.subfind));

% icm signifies that values are calculated over gas cells from main subhalo
% and fuzz, after excising other satellites. 
galFoFProps.icm=galFoFProps.full;

typeTag={'full','icm'};
%%run over simulations and snapshots
for k=1:length(sims)
    
    bp=illustris.set_env(sims(k));
    %global simDisplayName
    global LBox
    fprintf('running on simulation %s \n ', sims(k));
    
    simMask=strcmp(objectTable.sim,sims(k));
    
    for i=1:length(snaps)  % run over snaps 
        
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
        fprintf('looping over hosts \n');
        for j=1:length(hostList)
            
            hind=hostList(j);
            % find the relavant sat's in this particular host 
            localGalList=galInds(hostInds==hind);
            localTabInds=tabInds(hostInds==hind);  
            
            
            % load FOF gas for a host 
            gas=illustris.snapshot.loadHalo(bp,snap,hind-1,'gas',fields);
        
            % generate mask for ICM gas (main subhalo and 'fuzz')    
            icmMask= get_fof_icm(fofs,hind-1,snap);
            
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
            for jj=1:length(localGalList)
                
                % find index of galaxy in objectTable 
                locInd=localTabInds(jj);
                
                %center fof gas around galaxy
                locCoord=illustris.utils.centerObject(gas.CenterOfMass,subs.SubhaloPos(:,localGalList(jj)));
                gasDist=sqrt( sum(double(locCoord).^2,1)).*illUnits.lengthUnit;  %distance in physical kpc
                
                galVel=subs.SubhaloVel(:,localGalList(jj)).*illustris.utils.velocityFactor(snap,'sub');
                
                %% loop over aperatures
                for ii=1:length(aperatures)
                    
                    aper=apertures(ii);
                    for kk=1:length(typeTag)
                        
                    % generate mask around object
                    switch typeTag{kk}
                        case 'full'
                            mask=gasDist<=aper;
                        case 'icm'
                            mask=gasDist<=aper & icmMask;
                    end
                    
                    
                    
                    % local density 
                    rhoLocal=mean(gas.Density(mask));
                    galFoFProps.(typeTag{kk}).rhoLocal(ii,locInd)=rhoLocal;
                    
                    % local Velocity
                    for kk=1:3
                        velLocal(kk)=sum(gas.Masses(mask).*gas.newVel(kk,mask));
                    end
                    velLocal=velLocal./sum(gas.Masses(mask));
                    galFoFProps.(typeTag{kk}).velMedium(ii,:,locInd)=velLocal;
                    
                    % ram Pressure 
                    galFoFProps.(typeTag{kk}).ramPress(ii,locInd)=rhoLocal.*sum((galVel-velLocal).^2);
                    
                    % sonic velocity 
                    csonic=csFac.*sqrt(gas.Temperature(mask));
                    
                    
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
    
    
    
    
