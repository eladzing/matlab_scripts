%% get the scores and object list 

global DEFAULT_MATFILE_DIR
global simDisplayName
load([DEFAULT_MATFILE_DIR '/jf_objectTable_TNG50.mat'])
load([DEFAULT_MATFILE_DIR '/jellyfishScores_TNG50.mat'])


%% Generate catalogs by snapshots
 
snaps=unique(objectTable.snap);

% run over snap
outStruct=struct('done',[],'Score',[]);
for i=1:length(snaps)
    
    fprintf('running on snapshot %s \n ', num2str(snaps(i)));

    subs=illustris.groupcat.loadSubhalos(bp,snaps(i),{'SubhaloFlag'});
    
    galScore=zeros(length(subs),1)-1; % set all to -1 
    
    % find the objects from the right snapshot 
    snapMask=objectTable.snap==snaps(i); 
    tabInd=find(snapMask); % indices in objectTable from the right snpashot 
    
    % find their subfind index in the table 
    subfindInd=objectTable.subfind(snapMask)+1; % subfind is zero based, so index is increases by 1
    
    %prepare the scores 
    galScore(subfindInd)=tally6(tabInd); 
    outStruct.Score=galScore;
    outStruct.done=int8(galScore~=-1);    
    
    
    
    %% write to catalog 
    catName='jellyfish_flags';
        
        folder=['jellyfish_flags/' simDisplayName];
        
        illustris.utils.write_catalog(outStruct,snaps(i),'name',catName,...
            'path','default','folder',folder,'v');
    fprintf('===============\n');
end
    
    
    
    
    
    
    
    
    
    
