

bp=illustris.set_env(100);

list=[194947,197109,343639,371014,397867,453754,517859,524875];
treeType='SubLink_gal';

for i=1:length(galHistory)
    disp(i)
    tic
    treeFields={'SnapNum','SubfindID','SubhaloVel','SubhaloPos','SubhaloHalfmassRadType',...
    'SubhaloID','SubhaloMassInRadType','SubhaloSFRinRad','Group_M_Crit200','SubfindID',...
    'SubhaloHalfmassRadType','FirstSubhaloInFOFGroupID','GroupPos','Group_R_Crit200'};
tree(i) = illustris.sublink.loadTree(BASEPATH,99,list(i)-1,treeFields,true,treeType);
    toc
  
end


%% find redshifts 

global subBox

sbmin=subBox.center'-subBox.boxSize/2;
sbmax=subBox.center'+subBox.boxSize/2;
for i=1:length(galHistory)
    
    pos=tree(i).SubhaloPos;
    for k=1:3
        m(k,:)=pos(k,:)>sbmin(k) & pos(k,:)<sbmax(k);
    end
    mm=all(m,1);
    snaps=tree(i).SnapNum(mm);
    snapmax(i)=snaps(end);
     
    %zredMax(i)=zz(end);
    clear pos m mm zz snaps
end
zredMax=illustris.utils.snap2redshift(snapmax);
% global DRACOFLAG
% if DRACOFLAG
%     global DEFAULT_MATFILE_DIR
%     global simDisplayName
%     fname=sprintf('list_history_%s',simDisplayName);
%     %fnameBird=sprintf('list_history_birds_%s',simDisplayName);
%     
%     save([DEFAULT_MATFILE_DIR '/' fname],'galHistory','galBirdHistory','galProf','-v7.3')
%     fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
%     %fprintf(' *** Bird Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fnameBird]);
%     
% end
% fprintf(' *** DONE!  *** \n');
% 
