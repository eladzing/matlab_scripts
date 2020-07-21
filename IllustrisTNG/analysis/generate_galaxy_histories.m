% generate histories of centrals.

%global matFilePath
global simDisplayName
if ~skip2point
    
    sim='100';
    snap=99; %z=0
    zmax=3.0;
    
    bp=illustris.set_env(sim,'draco');
    
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
    
    
    massThresh=10^9.5; % threshold for *stellar* mass
    
    %% load subsinfo
    subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs,snap);
    units; % load general unit structure in cgs.
    
    %hostMass=fofs.Group_M_Crit200(subsInfo.hostFof+1).*illUnits.massUnit;
    massAllGals= subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit; % stellar mass within 2*rhalf
    
    massMask=massAllGals>massThresh;
    %hostMassMask=hostMass>hostMassThresh;
    % this mask selects all galaxies with dm component & stars, above *stellar* mass limit whose host has virials parameters
    
    galMask= subs.SubhaloFlag & subsInfo.isCentral & subsInfo.DM10perc & subsInfo.hasStars & massMask & subsInfo.hostHasVirial ;
    
    
    
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
id0=0;
nextP=0;

done=false(size(galMask));

for id=0:len-1
    if cnt>10
       break
    end
    
    perCent=floor((id+1)/len*100);
    
    if (perCent>=stepNext)     %       mod(perCent,5)==0 && perCent>=5)
        fprintf('%s %% of Galaxies done  \n',num2str(floor(perCent)))
        stepNext=stepNext+step;
    end
    
    
    if groupMask(id+1)
        cnt=cnt+1;
        
        [galHistory(id+1),~]=get_galaxy_history_full(id,'zmax',zmax,'startSnap',snap,'nobird');
        done(id+1)=true;
        
        
        prc=floor(cnt./sum(groupMask)*100);
        if (prc>=nextP)
            
            fprintf(' *** %i %% of histories completed  *** \n',prc);
            nextP=nextP+2;
            
        end
        
        %write interim
        if mod(prc,10) && prc>0
            fname=sprintf('galaxyHistory_z%s_%s_id%i',num2str(illustris.utils.get_zred(snap)),simDisplayName,id);
            save([DEFAULT_MATFILE_DIR '/' fname],'galHistory','done','-v7.3')
            fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
            
            if birdFlag
                fnameBird=sprintf('galaxyHistory_birds_z%s_%s_id%i',num2str(illustris.utils.get_zred(snap)),simDisplayName,id);
                save([DEFAULT_MATFILE_DIR '/' fnameBird],'galBirdHistory','done','-v7.3')
                fprintf(' *** Bird Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fnameBird]);
            end
        end
        
    end
end
if ~done(len)
    ii=find(~done,1,'first');
    galHistory(len)=galHistory(ii);
    galBirdHistory(len)=galBirdHistory(ii);
end


fname=sprintf('galaxyHistory_z%s_%s',num2str(illustris.utils.get_zred(snap)),simDisplayName);
save([DEFAULT_MATFILE_DIR '/' fname],'galHistory','done','-v7.3')
fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);

if birdFlag
    fnameBird=sprintf('galaxyHistory_birds_z%s_%s',num2str(illustris.utils.get_zred(snap)),simDisplayName);
    save([DEFAULT_MATFILE_DIR '/' fnameBird],'galBirdHistory','done','-v7.3')
    fprintf(' *** Bird Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fnameBird]);
end

%% organize output


% global DRACOFLAG
% if DRACOFLAG
%     global DEFAULT_MATFILE_DIR
%     global simDisplayName
%     fname=sprintf('BCG_history_z%s_%s',num2str(illustris.utils.get_zred(snap)),simDisplayName);
%     fnameBird=sprintf('BCGGroup_history_birds_z%s_%s',num2str(illustris.utils.get_zred(snap)),simDisplayName);
%     save([DEFAULT_MATFILE_DIR '/' fname],'galHistory','done','massThresh','-v7.3')
%     save([DEFAULT_MATFILE_DIR '/' fnameBird],'galBirdHistory','done','massThresh','-v7.3')
%     fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
%     fprintf(' *** Bird Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fnameBird]);
%
% end
% fprintf(' *** DONE!  *** \n');

