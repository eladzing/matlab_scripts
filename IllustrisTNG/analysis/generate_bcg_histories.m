% generate histories of centrals.

%global matFilePath
%global cosmoStruct
if ~skip2point
    
    sim='100';
    snap=99; %z=0
    zmax=4.0;
    
    bp=illustris.set_env(sim);
    
    global illUnits
    
    if ~exist('readFlag','var')
        readFlag=true;
    end
    
    if readFlag
        fprintf(' *** Reading data *** \n');
        fofs=illustris.groupcat.loadHalos(bp,snap);
        subs=illustris.groupcat.loadSubhalos(bp,snap);
        readFlag=false;
               
    end
    
    
    massThresh=10^9; % threshold for *stellar* mass
    hostMassThresh=10^(13.5);
    %% load subsinfo
    subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs,snap);
    units; % load general unit structure in cgs.
    
    hostMass=fofs.Group_M_Crit200(subsInfo.hostFof+1).*illUnits.massUnit;
    massAllGals= subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit; % stellar mass within 2*rhalf
       
    massMask=massAllGals>massThresh;
    hostMassMask=hostMass>hostMassThresh;
    % this mask selects all galaxies with dm component & stars, above *stellar* mass limit whose host has virials parameters
    galMask= subsInfo.isCentral & subsInfo.DM10perc & subsInfo.hasStars & massMask & hostMassMask & subsInfo.hostHasVirial & subsInfo.hasGas ;%  & subs.SubhaloMassInRadType(1,:)>0;
      
    
  
%     galMass=tCoolStruct.galMass;  % galaxy stellar mass
%     %gasMass=tCoolStruct.inGal.gasMass(1,:);
%     gasEnt=tCoolStruct.inGal.meanEntMW(1,:);
%     
%     
%     ssfr=illustris.utils.calc_ssfr(subs);
    
    %% select galaxies
    
   
    groupMask= galMask;
    %qGroupInd=find(groupMask);
    
end

%% generate values
fprintf(' *** Running over %s Galaxies *** \n',num2str(length(groupMask)));

step=5;
stepNext=5;
len=double(subs.count);
cnt=0;
nextP=0;

done=false(size(galMask));

for id=0:len-1
    
    perCent=floor((id+1)/len*100);
    
    if (perCent>=stepNext)     %       mod(perCent,5)==0 && perCent>=5)
        fprintf('%s %% of Galaxies done  \n',num2str(floor(perCent)))
        stepNext=stepNext+step;
    end
    
    
    if groupMask(id+1)
        cnt=cnt+1;
        
        prc=floor(cnt./sum(groupMask)*100);
        if (prc>=nextP)
            
            fprintf(' *** %i %% of histories completed  *** \n',prc);
            nextP=nextP+2;
        end
        
        [galHistory(id+1),galBirdHistory(id+1)]=get_galaxy_history_full(id,'zmax',zmax,'startSnap',snap);
        done(id+1)=true;
    end
end
if ~done(len)
    ii=find(~done,1,'first');
    galHistory(len)=galHistory(ii);
    galBirdHistory(len)=galBirdHistory(ii);
end

%% organize output


global DRACOFLAG
if DRACOFLAG
    global DEFAULT_MATFILE_DIR
    global simDisplayName
    fname=sprintf('BCG_history_z%s_%s',num2str(illustris.utils.get_zred(snap)),simDisplayName);
    fnameBird=sprintf('BCGGroup_history_birds_z%s_%s',num2str(illustris.utils.get_zred(snap)),simDisplayName);
    save([DEFAULT_MATFILE_DIR '/' fname],'galHistory','done','-v7.3')
    save([DEFAULT_MATFILE_DIR '/' fnameBird],'galBirdHistory','done','-v7.3')
    fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
    fprintf(' *** Bird Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fnameBird]);
    
end
fprintf(' *** DONE!  *** \n');

