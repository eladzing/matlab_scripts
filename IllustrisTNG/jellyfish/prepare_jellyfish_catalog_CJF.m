%% get the scores and object list

global DEFAULT_MATFILE_DIR

load([DEFAULT_MATFILE_DIR '/cosmic_jellyfish_objectTable.mat'])
%load([DEFAULT_MATFILE_DIR '/jellyfishScores_TNG50.mat'])

%% funky galaxies:
funkyTNG50.snap33=[21560, 50698, 90630, 127581, 138654];
funkyTNG50.snap40= [102, 126, 158, 36250, 66110, 112204, 163074];
funkyTNG50.snap50=79623;
funkyTNG50.snap59=[];
funkyTNG50.snap67=356638;
funkyTNG50.snap72=91931;
funkyTNG50.snap78=[226474, 451164];
funkyTNG50.snap84= [];
funkyTNG50.snap91= [];
funkyTNG50.snap99= [];

funkyTNG100.snap33=[];
funkyTNG100.snap40= [];
funkyTNG100.snap50=10011;
funkyTNG100.snap59=[301337, 301338];
funkyTNG100.snap67= [];
funkyTNG100.snap72= [];
funkyTNG100.snap78= [];
funkyTNG100.snap84= [];
funkyTNG100.snap91= [];
funkyTNG100.snap99= [];
%% Generate catalogs by simulations and snapshots snapshots based on an object Table


snaps=unique(objectTable.snap);
sims=unique(objectTable.sim);

% run over simulations
for k=1:length(sims)
    bp=illustris.set_env(sims(k),'nomount');
    global simDisplayName
    fprintf('running on simulation %s \n ', sims(k));
    
    simMask=strcmp(objectTable.sim,sims(k));
    
    switch sims(k)
        case 'TNG100'
            funky=funkyTNG100;
        case 'TNG50'
            funky=funkyTNG50;
    end
            
    
    % run over snap
    for i=1:length(snaps)
        outStruct=struct('done',[],'Score',[],'ClassificationNum',[],'ScoreTotal',[]);
        
        
        fprintf('running on snapshot %s \n ', num2str(snaps(i)));
        
        subs=illustris.groupcat.loadSubhalos(bp,snaps(i),{'SubhaloFlag'});
        
        galScore=zeros(length(subs),1)-1; % set all to -1
        galNcls=zeros(length(subs),1)-1;
        galScoreTot=zeros(length(subs),1)-1;
        
        % find the objects from the right snapshot
        snapMask=objectTable.snap==snaps(i);
        fullMask=snapMask & simMask;
        tabInd=find(fullMask); % indices in objectTable from the right simulation and snpashot
        
        % find their subfind index in the table
        subfindInd=objectTable.subfind(fullMask)+1; % subfind is zero based, so index is increases by 1
        
        %prepare the scores
        galScore(subfindInd)=objectTable.score(tabInd);
        galNcls(subfindInd)=objectTable.clsNum(tabInd);
        galScoreTot(subfindInd)=objectTable.scoreTotal(tabInd);
        
        %% address funky galaxyis 
        funkyID=funky.(['snap' num2str(snaps(i))]);
        
        if ~isempty(funkyID)
            galScore(funkyID+1)=0;
            galNcls(funkyID+1)=0;
            galScoreTot(funkyID+1)=0;
        end
        
        outStruct.Score=galScore;
        outStruct.ClassificationNum=galNcls;
        outStruct.ScoreTotal=galScoreTot;
        outStruct.done=int8(galScore~=-1);
        
        
        
        %% write to catalog
        catName='cosmic_jellyfish_flags';
        
        folder=['jellyfish_flags/' simDisplayName];
        
        illustris.utils.write_catalog(outStruct,snaps(i),'name',catName,...
            'path','default','folder',folder,'v');
        fprintf('===============\n');
    end
    
end  
    
    
    
    
    
    
    
    
