% generate histories of centrals.

%global matFilePath
%global cosmoStruct

if ~skip2point
    
    %sim='300';
    snap=99; %z=0
    zmax=6.0;
    
    %bp=illustris.set_env(sim);
    
    global illUnits
    global BASEPATH
    global LBox
    global simDisplayName
    %treeType='SubLink_gal';
    
    
    if ~exist('readFlag','var')
        readFlag=true;
    end
    
    
    if readFlag
        fprintf(' *** Reading data *** \n');
        fofs=illustris.groupcat.loadHalos(bp,snap);
        subs=illustris.groupcat.loadSubhalos(bp,snap);
        readFlag=false;
               
    end
    
    % assign a lower stellar mass threshold based on the TNG box. 
    if contains(simDisplayName,'100') || contains(simDisplayName,'300')
        massThresh=10^9; % threshold for *stellar* mass
    elseif contains(simDisplayName,'50')
        massThresh=10^7; % threshold for *stellar* mass
    else
        error('mass threshold not set  - could not identify simulation: %s \n',simDisplayName);
    end
    fprintf('Setting stellar mass threshold to: %0.1e solar mass \n',massThresh); 
    
    %hostMassThresh=10^(13.5);
    
    %% load subsinfo
    subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
    units; % load general unit structure in cgs.
    
    hostMass=fofs.Group_M_Crit200(subsInfo.hostFof+1).*illUnits.massUnit;
    massAllGals= subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit; % stellar mass within 2*rhalf

    % Generate mask for centrals from within the 'good' galaxies' 
    galMask=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'centrals');

    treeType='SubLink_gal';
    snapStart=snap;
    
    
end

%% run over galaxies and generate the histories of paramters 
fprintf(' *** Running over %s Galaxies *** \n',num2str(length(galMask)));

step=5;
stepNext=5;
len=double(subs.count);
cnt=0;
nextP=0;

done=false(size(galMask));
%centralHist=[];
for id=0:len-1
    
    perCent=floor((id+1)/len*100);
    
    if (perCent>=stepNext)     %       mod(perCent,5)==0 && perCent>=5)
        fprintf('%s %% of Galaxies done  \n',num2str(floor(perCent)))
        stepNext=stepNext+step;
    end
    
    
    if galMask(id+1)
        cnt=cnt+1;
        
        prc=floor(cnt./sum(galMask)*100);
        if (prc>=nextP)
            
            fprintf(' *** %i %% of histories completed  *** \n',prc);
            nextP=nextP+2;
        end
        
        
        
        % read trees
        treeFields={'SnapNum','SubfindID','SubhaloPos','SubhaloSFRinRad',...
            'SubhaloID','SubhaloMassInRadType','Group_M_Crit200',...
            'GroupFirstSub','GroupPos','Group_R_Crit200'};
        tree = illustris.sublink.loadTree(BASEPATH,snapStart,id,treeFields,true,treeType);
        
        redshifts=illustris.utils.snap2redshift(tree.SnapNum);
        [~,lastInd]=min(abs(redshifts-zmax));
               
        unitFactors=illustris.utils.calc_illUnits(tree.SnapNum(1:lastInd));

        %% Get things which can be gotten directly from tree. 
        res.hostMass=tree.Group_M_Crit200(1:lastInd)'.*unitFactors.massUnit;
        res.hostR200=tree.Group_R_Crit200(1:lastInd)'.*unitFactors.lengthUnit;
        %res.halfMassRadType=tree.SubhaloHalfmassRadType(1:lastInd).*unitFactors.lengthUnit;
        
        res.SubfindID=tree.SubfindID(1:lastInd)';
        
        res.sfr=tree.SubhaloSFRinRad(1:lastInd)';
        res.stellarMass=tree.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,1:lastInd).*unitFactors.massUnit;
        res.bhMassInRad=tree.SubhaloMassInRadType(illustris.partTypeNum('bh')+1,1:lastInd).*unitFactors.massUnit;
%         res.ssfr=res.sfr./res.stellarMass+...
%             10^0.5*1e-17.*10.^(0.5.*rand(size(res.sfr)));
        
        %res.isCentral=tree.SubhaloID(1:lastInd)==tree.FirstSubhaloInFOFGroupID(1:lastInd);
        res.isCentral=tree.SubfindID(1:lastInd)==tree.GroupFirstSub(1:lastInd);
        res.radiusToHost=findDistance(tree.SubhaloPos(:,1:lastInd),tree.GroupPos(:,1:lastInd),LBox,3).*unitFactors.lengthUnit;
        %aexp=illustris.utils.snap2aexpn(tree.SnapNum(1:lastInd));
        res.zred=illustris.utils.snap2redshift(tree.SnapNum(1:lastInd));
        %res.zredExt=res.zred(fullSnapMask(tree.SnapNum(1:lastInd)+1));
        
        %[galHistory(id+1),galBirdHistory(id+1)]=get_galaxy_history_full(id,'zmax',zmax,'startSnap',snap);
        centralHist(id+1)=res;
        done(id+1)=true;
    end
end
if ~done(len)
    ii=find(~done,1,'first');
    centralHist(len)=centralHist(ii);
    
end

%% organize output


global DRACOFLAG
if DRACOFLAG
    global DEFAULT_MATFILE_DIR
    
    fname=sprintf('central_history_splashback_snp%s_%s',num2str(snap),simDisplayName);
    
    save([DEFAULT_MATFILE_DIR '/' fname],'centralHist','done','-v7.3')
    
    fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
    
    
end
fprintf(' *** DONE!  *** \n');

