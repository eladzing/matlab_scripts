%% track a subhalo in the sub-box based on some initial timestep
tic
sim='100';
bp=illustris.set_env(sim);
%global BASEPATH;
global DEFAULT_MATFILE_DIR

zEnd=4;


load([DEFAULT_MATFILE_DIR '/sub0_cgmList_TNG100.mat'])

treeId=1;
galID=tree(treeId).SubfindID(1); % 194946;

snaps=tree(treeId).SnapNum;
zredTNG=illustris.utils.snap2redshift(snaps);


%zred0=illustris.utils.snap2redshift(startSnapTNG:-1:endSnapTNG);

%% sub box stuff
bpSub=illustris.set_env(sim,'sub0');
global subBox
snapSub=subBox.Nsnaps-1:-1:0;
zredSub=illustris.utils.snap2redshift(snapSub);

% find when subhalo leaves the subbox
sbmin=subBox.center'-subBox.boxSize/2;
sbmax=subBox.center'+subBox.boxSize/2;
sbMask=tree(treeId).SubhaloPos(1,:)>sbmin(1) & tree(treeId).SubhaloPos(1,:)<sbmax(1) & ...
    tree(treeId).SubhaloPos(2,:)>sbmin(2) & tree(treeId).SubhaloPos(2,:)<sbmax(2) & ...
    tree(treeId).SubhaloPos(3,:)>sbmin(3) & tree(treeId).SubhaloPos(3,:)<sbmax(3);

if any(~sbMask)
    occInd=find(~sbMask,1,'first');
    zEnd=min(zEnd,zredTNG(occInd));
end

% find final redshift
zIndTNG=1:find(zredTNG>zEnd,1,'first');
zIndSub=1:find(zredSub>zredTNG(zIndTNG(end)),1,'first')-1;

% identify when bigbox snapshots are found
mask=true(size(zIndSub));
cnt=1;
for i=zIndSub
    if zredTNG(cnt)==zredSub(i)
        mask(i)=false;
        cnt=cnt+1;
    end
end

%% run back and get center

subCenter=zeros(3,zIndSub(end));
fac=0.05;

cnt=0;
step1=1;
step=1;
skipCount=0;
skip=20;
for ii=zIndSub
    prc=100.*ii./length(zIndSub);
    if prc>step
        fprintf(' %s %% finished \n',num2str(floor(prc)))
        step=step+step1;
    end
    
    if ~mask(ii)    % big box snapshot exists
        cnt=cnt+1;
        zi=zIndTNG(cnt);
        subCenter(:,ii)=tree(treeId).SubhaloPos(:,zi);
        
        %read stars
        starIds=illustris.snapshot.loadSubhalo(bp,snaps(zi),tree(treeId).SubfindID(zi),illustris.partTypeNum('stars'),{'ParticleIDs'});
        
        
        len=ceil(fac*length(starIds));
        starIds=starIds(1:len);
        
        
        
    else
        skipCount=skipCount+1;
        
        if skipCount==skip
            skipCount=0;
            
            stars=illustris.snapshot.loadSubset(bpSub,snapSub(ii),illustris.partTypeNum('stars'),{'Coordinates','Masses','ParticleIDs'});
            
            
            %find the sub stars identical to subfind stars
            matchedInd=zeros(1,len);
            for i=1:length(starIds)
                kk=find(stars.ParticleIDs==starIds(i));
                if ~isempty(kk)
                    matchedInd(i)=kk;
                end
            end
            
            
            % find center
            matchedInd=matchedInd(matchedInd~=0);
            
            for k=1:3
                subCenter(k,ii)=sum(stars.Coordinates(k,matchedInd).*stars.Masses(matchedInd))./sum(stars.Masses(matchedInd));
            end
            
        end
        
    end
end
toc

%% save

illustris.set_env(sim,'sub0');
subSnap=snapSub(zIndSub);
zFinal=illustris.utils.snap2redshift(subSnap(end));
%global DEFAULT_MATFILE_DIR
global simDisplayName
fname=sprintf('cgmList_sub%s_posHistory_z%1.2f_%s',num2str(galID),zFinal,simDisplayName);
save([DEFAULT_MATFILE_DIR '/' fname],'subSnap','subCenter','-v7.3')
fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);




% %% go back in time
%
% illustris.set_env(100,'sub0');
% global subBox
% startSnapSub=subBox.Nsnaps-1;
%
%
% zred=illustris.utils.snap2redshift(0:subBox.Nsnaps-1);
%
% %% match redshifts
%
% zInd=zeros(size(zred0));
%
% for i=1:length(zInd)
%     dz=abs(zred-zred0(i));
%     [~,zInd(i)]=min(dz);
% end
%
% %% match stars
% fac=0.05;
% len=ceil(fac*stars0.count);
% ids=stars0.ParticleIDs(1:len);
%
%
%
% cnt=0
% for j=zInd(1):-1:zInd(2)
%
%     cnt=cnt+1;
%     stars=illustris.snapshot.loadSubset(bp,j-1,illustris.partTypeNum('stars'));
%
%
% ind=zeros(1,len);
%
% %find the sub stars identical to subfind stars
% for i=1:length(ids)
%     ii=find(stars.ParticleIDs==ids(i));
%     if ~isempty(ii)
%         ind(i)=ii;
%     end
% end
%
%
%
% % find center
% ind=ind(ind~=0);
% for k=1:3
%     com(k,i)=sum(stars.Coordinates(k,ind).*stars.Masses(ind))./sum(stars.Masses(ind));
% end
%
%
%