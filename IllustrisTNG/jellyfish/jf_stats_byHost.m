

global illUnits
global DEFAULT_MATFILE_DIR

load([DEFAULT_MATFILE_DIR '/cosmic_jellyfish_objectTable.mat'])
load([DEFAULT_MATFILE_DIR '/jf_galProperties_CJF.mat']);

snaps=unique(objectTable.snap);
sims=unique(objectTable.sim);


%% run over simulations and snapshots
indx=0; %running index of hosts
for k=1:length(sims)
    
    bp=illustris.set_env(sims(k));
    %global simDisplayName
    %     global LBox
    fprintf('running on simulation %s \n ', sims(k));
    
    fprintf('going over %i snapshots \n',length(snaps));
    simMask=strcmp(objectTable.sim,sims(k));
    
    for i=1:length(snaps)
        
        snap=snaps(i);
        fprintf('%i )  snap %i ',i,snap);
        % find indices of gals from the snapshot
        snapMask=objectTable.snap==snap;
        %fullMask=snapMask & simMask;
        
        inds=find(snapMask & simMask); % indices of relavant object in the table
        hostList=unique(galProps.hostID(inds)); %id's (0-based) of hosts
        hostInds=hostList+1; % indices in catalog of hosts
        
        if isempty(inds) % skip snapshot with no relevant objects
            fprintf(' .... no objects \n');
            continue
        else
            fprintf(' ... %i hosts \n',length(hostList));
        end
        
        fprintf('   Loading catalogs of snap %i \n',snap);
        [subs,fofs,subsInfo]=illustris.loadFofSub(snap);
        illustris.utils.set_illUnits(snap);
        
        fprintf('   running over hosts \n');
        
        
        
        for j=1:length(hostList)
            indx=indx+1;
            
            tag=join([sims(k) 'snp' num2str(snap,'%03.f') 'hostID' num2str(hostList(j))],'');
            
            jfStats.GroupNsubs(indx)=fofs.GroupNsubs(hostInds(j));
            jfStats.sim(indx)=sims(k);
            jfStats.tag(indx)=tag;
            jfStats.snap(indx)=snap;
            jfStats.M200c(indx)=double(fofs.Group_M_Crit200(hostInds(j))).*illUnits.massUnit;
            jfStats.hostID(indx)=hostList(j);
            
                        % find the relavent subs
            satList=find(snapMask & simMask & galProps.hostID'==hostList(j)); %indices of subs in the objectTable
            
            jfStats.sampleSats(indx)=length(satList);
            jfStats.JFNumRaw(indx)=sum(objectTable.scoreRaw(satList)>15);
            jfStats.JFNumWeighted(indx)=sum(objectTable.scoreWeighted(satList)>0.8);
            
            % fill out host info
            jfStats.halo(indx).tag=tag;
            jfStats.halo(indx).satTags=objectTable.tag(satList);
            jfStats.halo(indx).subfind=objectTable.subfind(satList);
            jfStats.halo(indx).satIndxInTable=satList;
            jfStats.halo(indx).scoreRaw=objectTable.scoreRaw(satList);
            jfStats.halo(indx).scoreWeighted=objectTable.scoreWeighted(satList);
            %jfStats.halo(indx).rpos=galProps.rpos(satList)./galProps.hostR200c(satList);
            
            
        end
        
        
        %         galInds=objectTable.subfind(inds)+1; % indices in the catalogs
        %         hostInds=subsInfo.hostFof(galInds)+1;% indices in the catalogs
        
        
        
    end
end

fname=sprintf('jf_statsByHosts_CJF.mat');
save([DEFAULT_MATFILE_DIR '/' fname],'jfStats','-v7.3')

fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);


%      %% fill out the properties
%
%
%
%         % virial properties
%         jfStats.hostID(inds)=subsInfo.hostFof(galInds);
%
%         jfStats.hostM200c(inds)=illUnits.massUnit.*...
%             double(fofs.Group_M_Crit200(hostInds));
%
%         jfStats.hostR200c(inds)=illUnits.lengthUnit.*...
%             double(fofs.Group_R_Crit200(hostInds));
%
%         jfStats.hostM200m(inds)=illUnits.massUnit.*...
%             double(fofs.Group_M_Mean200(hostInds));
%
%         jfStats.hostR200m(inds)=illUnits.lengthUnit.*...
%             double(fofs.Group_R_Mean200(hostInds));
%
%         % gal/subfind properties
%         jfStats.stellarMass(inds)=illUnits.massUnit.*...
%             double(subs.SubhaloMassType(illustris.partTypeNum('stars')+1,galInds));
%
%         jfStats.galStellarMass(inds)=illUnits.massUnit.*...
%             double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,galInds));
%
%         jfStats.gasMass(inds)=illUnits.massUnit.*...
%             double(subs.SubhaloMassType(illustris.partTypeNum('gas')+1,galInds));
%
%         jfStats.galGasMass(inds)=illUnits.massUnit.*...
%             double(subs.SubhaloMassInRadType(illustris.partTypeNum('gas')+1,galInds));
%
%         jfStats.galBHMass(inds)=illUnits.massUnit.*...
%             double(subs.SubhaloMassInRadType(illustris.partTypeNum('BH')+1,galInds));
%
%         jfStats.galSFR(inds)=...
%             double(subs.SubhaloSFRinRad(galInds));
%
%         jfStats.SFR(inds)=...
%             double(subs.SubhaloSFR(galInds));
%
%         %% position data
%         galPos=double(subs.SubhaloPos(:,galInds)); % global position in simulation box, in simulation units
%         hostPos=double(fofs.GroupPos(:,hostInds));% global position of host in simulation box, in simulation units
%
%
%         % apply periodic boundary conditions to calculate the actual radial
%         % position
%         jfStats.pos(:,inds)=galPos;
%
%         rpos=findDistance(galPos,hostPos,LBox,3);
%         jfStats.rpos(inds)=rpos.*illUnits.lengthUnit;
%
%
%
%         % leave this for now till we test.
%         for ii=1:3
%             mask=abs(galPos(ii,:)-hostPos(ii,:))>0.5.*LBox;
%
%             if any(mask)
%                 mask1=mask & hostPos(ii,:)>0.5*LBox;
%                 mask2=mask & hostPos(ii,:)<=0.5*LBox;
%
%                 galPos(ii,mask1)=galPos(ii,mask1)+LBox;
%                 galPos(ii,mask2)=galPos(ii,mask2)-LBox;
%             end
%
%             galPos(ii,:)=galPos(ii,:)-hostPos(ii,:);
%
%         end
%
%         % velocity w.r.t host
%         vsat=double(subs.SubhaloVel(:,galInds)).*illustris.utils.velocityFactor(snap,'sub'); % 3d velocity for each galaxy
%         vhost=double(fofs.GroupVel(:,hostInds)).*illustris.utils.velocityFactor(snap,'host');% 3d velocity of host for each galaxy
%
%         vel=vsat-vhost;
%
%
%
%         %        vrad=sum(vel.*galPos,1)./rpos;
%
%         jfStats.vel(:,inds)=vel;
%         jfStats.vrad(inds)=sum(vel.*galPos,1)./rpos;
%
%
