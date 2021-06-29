% generate histories of centrals.

%global matFilePath
%global cosmoStruct
hostID=0;
if ~skip2point
    
    %sim='100';
    snap=99; %z=0
    zmax=4.0;
    
    %bp=illustris.setgen_env(sim);
    
    global illUnits
    
    if ~exist('readFlag','var')
        readFlag=true;
    end
    
    if readFlag
        fprintf(' *** Reading data *** \n');
        loadFofSub
        
        readFlag=false;
        
    end
    
    massThresh=10^9; % threshold for *stellar* mass
    hostMassThresh=10^(13.5);
    %% load subsinfo
    subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
    units; % load general unit structure in cgs.
    
    hostMass=fofs.Group_M_Crit200(subsInfo.hostFof+1).*illUnits.massUnit;
    massAllGals= illustris.utils.get_stellar_mass(subs);%  subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit; % stellar mass within 2*rhalf
    
    centralMask=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'centrals');
    satMask=    illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'sats');
    
    hostMask=subsInfo.hostFof==hostID;%     hostMass>hostMassThresh;
    
    ssfr=illustris.utils.calc_ssfr(subs);
    
    %% select galaxies
    
    
    groupMask= satMask & hostMask;
    %qGroupInd=find(groupMask);
    
end

%% generate values
fprintf(' *** Running over %s Galaxies *** \n',num2str(sum(groupMask)));

satID=find(groupMask)-1;

step=5;
stepNext=5;
len=length(satID);
cnt=0;
nextP=0;

%done=false(size(galMask));

for k=1:length(satID)
    
    %     perCent=floor((id+1)/len*100);
    %
    %     if (perCent>=stepNext)     %       mod(perCent,5)==0 && perCent>=5)
    %         fprintf('%s %% of Galaxies done  \n',num2str(floor(perCent)))
    %         stepNext=stepNext+step;
    %     end
    %
    
    %if groupMask(id+1)
    %   cnt=cnt+1;
    
    prc=floor(k./sum(groupMask)*100);
    if (prc>=nextP)
        
        fprintf(' *** %i %% of histories completed  *** \n',prc);
        nextP=nextP+2;
    end
    
    [galHistory(k),galBirdHistory(k)]=get_galaxy_history_full(satID(k),'zmax',zmax,'startSnap',snap);
    %done(id+1)=true;
end

% if ~done(len)
%     ii=find(~done,1,'first');
%     galHistory(len)=galHistory(ii);
%     galBirdHistory(len)=galBirdHistory(ii);
% end

%% organize output


global DRACOFLAG
if DRACOFLAG
    global DEFAULT_MATFILE_DIR
    global simDisplayName
    fname=sprintf('satFof0_history_0%s_%s',num2str(snap),simDisplayName);
    fnameBird=sprintf('satFof0_birdHhistory_0%s_%s',num2str(snap),simDisplayName);
    save([DEFAULT_MATFILE_DIR '/' fname],'galHistory','-v7.3')
    save([DEFAULT_MATFILE_DIR '/' fnameBird],'galBirdHistory','-v7.3')
    fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
    fprintf(' *** Bird Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fnameBird]);
    
end
fprintf(' *** DONE!  *** \n');

