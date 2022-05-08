

global illUnits
global DEFAULT_MATFILE_DIR

load([DEFAULT_MATFILE_DIR '/cosmic_jellyfish_objectTable.mat'])
snaps=unique(objectTable.snap);
sims=unique(objectTable.sim);
%% initialize

% snap & sim and id
galProps.snap=objectTable.snap;
galProps.sim=objectTable.sim;
galProps.subfind=objectTable.subfind;
galProps.tag=objectTable.tag;

ngal=length(objectTable.subfind);



% galProps.sim=zeros(1,ngal);
% galProps.snap=zeros(1,ngal);
galProps.stellarMass=zeros(1,ngal);
galProps.hostM200c=zeros(1,ngal);
galProps.hostR200c=zeros(1,ngal);
galProps.hostM200m=zeros(1,ngal);
galProps.hostR200m=zeros(1,ngal);
galProps.galStellarMass=zeros(1,ngal);
galProps.gasMass=zeros(1,ngal);  % total
galProps.galGasMass=zeros(1,ngal); % total 
galProps.galSFR=zeros(1,ngal);
galProps.pos=zeros(3,ngal);
galProps.vel=zeros(3,ngal);
galProps.rpos=zeros(1,ngal);
galProps.vrad=zeros(1,ngal);
%galProps.vtan=zeros(1,ngal);
galProps.galBHMass=zeros(1,ngal);
galProps.hostID=zeros(1,ngal);
%galProps.nSatHost=zeros(1,ngal);
%galProps.rhoLocal=zeros(size(objectTable.subfind));
%galProps.ramPress=zeros(size(objectTable.subfind));
%galProps.mach=zeros(size(objectTable.subfind));

%%run over simulations and snapshots
for k=1:length(sims)
    
    bp=illustris.set_env(sims(k));
    %global simDisplayName
    global LBox
    fprintf('running on simulation %s \n ', sims(k));
    
    fprintf('going over %i snapshots \n',length(snaps));
    simMask=strcmp(objectTable.sim,sims(k));
    
    for i=1:length(snaps)
        
        snap=snaps(i);
        fprintf('%i )  snap %i ',i,snap);
        % find indices of gals from the snapshot
        snapMask=objectTable.snap==snap;
        %fullMask=snapMask & simMask;
        inds=find(snapMask & simMask);
        
        if isempty(inds) % skip snapshot with no relevant objects
            fprintf(' .... no objects \n');
            continue
        else
            fprintf(' ... %i objects \n',length(inds));
        end
        
        fprintf('   Loading catalogs of snap %i \n',snap);
        [subs,fofs,subsInfo]=illustris.loadFofSub(snap);
        illustris.utils.set_illUnits(snap);
                
        fprintf('   Getting galaxy properties \n');
        
        
        
        galInds=objectTable.subfind(inds)+1; % indices in the catalogs 
        hostInds=subsInfo.hostFof(galInds)+1;% indices in the catalogs 
        
        %% fill out the properties
        
        
        
        % virial properties
        galProps.hostID(inds)=subsInfo.hostFof(galInds);
        
        galProps.hostM200c(inds)=illUnits.massUnit.*...
             double(fofs.Group_M_Crit200(hostInds));
        
        galProps.hostR200c(inds)=illUnits.lengthUnit.*...
             double(fofs.Group_R_Crit200(hostInds));
        
        galProps.hostM200m(inds)=illUnits.massUnit.*...
             double(fofs.Group_M_Mean200(hostInds));
        
        galProps.hostR200m(inds)=illUnits.lengthUnit.*...
             double(fofs.Group_R_Mean200(hostInds));
         
        % gal/subfind properties
        galProps.stellarMass(inds)=illUnits.massUnit.*...
             double(subs.SubhaloMassType(illustris.partTypeNum('stars')+1,galInds));
        
        galProps.galStellarMass(inds)=illUnits.massUnit.*...
             double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,galInds));
        
        galProps.gasMass(inds)=illUnits.massUnit.*...
             double(subs.SubhaloMassType(illustris.partTypeNum('gas')+1,galInds));
        
        galProps.galGasMass(inds)=illUnits.massUnit.*...
             double(subs.SubhaloMassInRadType(illustris.partTypeNum('gas')+1,galInds));
        
        galProps.galBHMass(inds)=illUnits.massUnit.*...
             double(subs.SubhaloMassInRadType(illustris.partTypeNum('BH')+1,galInds));
        
        galProps.galSFR(inds)=...
             double(subs.SubhaloSFRinRad(galInds));
         
         galProps.SFR(inds)=...
             double(subs.SubhaloSFR(galInds));
        
       %% position data 
       galPos=double(subs.SubhaloPos(:,galInds)); % global position in simulation box, in simulation units 
       hostPos=double(fofs.GroupPos(:,hostInds));% global position of host in simulation box, in simulation units 
       
       
       % apply periodic boundary conditions to calculate the actual radial
       % position 
       galProps.pos(:,inds)=galPos; 
       
       rpos=findDistance(galPos,hostPos,LBox,3);
       galProps.rpos(inds)=rpos.*illUnits.lengthUnit;
       
       
       
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
      
       galProps.vel(:,inds)=vel;
       galProps.vrad(inds)=sum(vel.*galPos,1)./rpos;
       
        
    end
end

%% add a uniqaue halo tag identifier 
galProps = jellyfish.generate_hostTag(galProps);



%% write to file 
fname=sprintf('jf_galProperties_CJF.mat');
save([DEFAULT_MATFILE_DIR '/' fname],'galProps','-v7.3')

fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);



    
    
    % galProps.vel(inds,3)=...
    %
    %     galProps.rpos(inds)=...
    %
    %     galProps.vrad(inds)=...
    %
    %     galProps.vtan(inds)=...
    
    %massThresh=10^8.5; % threshold for *stellar* mass
    
    
    
    
