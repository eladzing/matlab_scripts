%% Need build an objectTable similar to the CJF one but which includes all
% satellites which fit the CJF project criteria

snap100=[33 40 50 59 67 72 78 84 91 99];
snap50= [33 40 50 59 67:99];

global illUnits
global DEFAULT_MATFILE_DIR

load([DEFAULT_MATFILE_DIR '/cosmic_jellyfish_objectTable.mat'])
%snaps=unique(objectTable.snap);
sims=unique(objectTable.sim);

%% Initialize  table

% vnames={'sim','snap','subfind','hostID',...
%     'hostM200c','hostR200c','hostM200m','hostR200m',...
%     'stellarMass','galStellarMass','gasMass','galGasMass',...
%     'galSFR','galBHMass','rpos','vrad'};
% %varTypes=vnames;
% %varTypes{1}='string';
% varTypes = {'string','double','double','double','double','double',...
%     'double','double','double','double','double',...
%     'double','double','double','double','double'};
% fullSat=table([],[],[],[],[],[],[],[],[],[],...
%     [],[],[],[],[],[],'VariableTypes',varTypes,...
%     'variableNames',vnames);
% 
indStart=1;
for k=1:length(sims)
    
    %bp=illustris.set_env(sims(k));
    bp=illustris.set_env(sims(k),'nomount');
    
    if sims(k)=="TNG100"
        massThresh=10^9.5;
        snaps=snap100;
    elseif sims(k)=="TNG50"
        massThresh=10^8.3;
        snaps=snap50;
    else
        error('%s - unknown simulation: %s',current_function().upper,sims(k));
    end
    global LBox
    
    
    for i=length(snaps)
        %for i=1:length(snaps)
        snap=snaps(i);
        illustris.utils.set_illUnits(snap);
        
        %% load catalogs
        [subs,fofs,subsInfo]=illustris.loadFofSub(snap);
        
        %% generate mask
        smass=illustris.utils.get_stellar_mass(subs);
        
        msk=smass>=massThresh & ~subsInfo.isCentral & subs.SubhaloFlag;
        
        % check that the mask encompass the CJF sample
        sampleMask=objectTable.sim==sims(k) & objectTable.snap==snap;
        sampleID=objectTable.subfind(sampleMask);
        maskID=find(msk)-1;
        missingInds=[];
        for j=1:length(sampleID)
            if ~any(sampleID(j)==maskID)
                warning('%s - missing index from sample(%s, %i) : %i',...
                    current_function().upper,sims(k),snap,sampleID(j));
                missingInds(end+1)=sampleID(j);
            end
        end
        
        %% fill up table
        indEnd=indStart+sum(msk)-1;
        tabInds=indStart:indEnd;
        
        galInds=maskID+1; % indices in the catalogs
        hostInds=subsInfo.hostFof(galInds)+1;
        
        
        fullSat.sim(tabInds)=sims(k);
        fullSat.snap(tabInds)=snap;
        fullSat.subfind(tabInds)=maskID;
        fullSat.hostID(tabInds)=subsInfo.hostFof(msk);
        
        %     % build tag
        %     snp=string(['0' num2str(snap)]);
        %     tg1(1:length(mskInds),1)="snp";
        %     tg2(1:length(mskInds),1)="subid";
        %     tg3(1:length(mskInds),1)="typ:";
        %     clsTab.tag=join(['TNG100' tg1 snp tg2 sid tg3 imgType'],'');
        %     fullSat.tag(tabInds)=
        
        %% fill up properties
        
        
        
        fullSat.hostM200c(tabInds)=illUnits.massUnit.*...
            double(fofs.Group_M_Crit200(hostInds));
        fullSat.hostR200c(tabInds)=illUnits.lengthUnit.*...
            double(fofs.Group_R_Crit200(hostInds));
        fullSat.hostM200m(tabInds)=illUnits.massUnit.*...
            double(fofs.Group_M_Mean200(hostInds));
        fullSat.hostR200m(tabInds)=illUnits.lengthUnit.*...
            double(fofs.Group_R_Mean200(hostInds));
        fullSat.stellarMass(tabInds)=illUnits.massUnit.*...
            double(subs.SubhaloMassType(illustris.partTypeNum('stars')+1,galInds));
        fullSat.galStellarMass(tabInds)=illUnits.massUnit.*...
            double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,galInds));
        fullSat.gasMass(tabInds)=illUnits.massUnit.*...
            double(subs.SubhaloMassType(illustris.partTypeNum('gas')+1,galInds));
        fullSat.galGasMass(tabInds)=illUnits.massUnit.*...
            double(subs.SubhaloMassInRadType(illustris.partTypeNum('gas')+1,galInds));
        fullSat.galSFR(tabInds)=double(subs.SubhaloSFRinRad(galInds));
        fullSat.galBHMass(tabInds)=illUnits.massUnit.*...
            double(subs.SubhaloMassInRadType(illustris.partTypeNum('BH')+1,galInds));
        
        %% position and velocity data
        galPos=double(subs.SubhaloPos(:,galInds)); % global position in simulation box, in simulation units
        hostPos=double(fofs.GroupPos(:,hostInds));% global position of host in simulation box, in simulation units
        rpos=findDistance(galPos,hostPos,LBox,3);
        
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
        %fullSat.pos=
        %fullSat.vel(tabInds)=vel;
        fullSat.rpos(tabInds)=rpos.*illUnits.lengthUnit;
        fullSat.vrad(tabInds)=sum(vel.*galPos,1)./rpos;
        
        indStart=indEnd+1;
    end
end

%% write to file 

% vnames={'sim','snap','subfind','hostID',...
%     'hostM200c','hostR200c','hostM200m','hostR200m',...
%     'stellarMass','galStellarMass','gasMass','galGasMass',...
%     'galSFR','galBHMass','rpos','vrad'};
% %varTypes=vnames;
% varTypes = {'string','double','double','int32','double','double',...
%      'double','double','double','double','double',...
%      'double','double','double','double','double'};
% fullSatTable=table(fullSat.sim',fullSat.snap',fullSat.subfind',fullSat.hostID',...
%     fullSat.hostM200c',fullSat.hostR200c',fullSat.hostM200m',fullSat.hostR200m',...
%     fullSat.stellarMass',fullSat.galStellarMass',fullSat.gasMass',fullSat.galGasMass',...
%     fullSat.galSFR',fullSat.galBHMass',fullSat.rpos',fullSat.vrad',...
%     'VariableTypes',varTypes,'variableNames',vnames);
% 



fname=sprintf('jf_fullSatSample_CJF.mat');
save([DEFAULT_MATFILE_DIR '/' fname],'fullSat','-v7.3')

fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);



    