%% get the scores and object list

global DEFAULT_MATFILE_DIR

load([DEFAULT_MATFILE_DIR '/cosmic_jellyfish_objectTable.mat'])


%% Generate a single catalog containing all objects with a random and preferred 
% image orientation. The catalog must contain sim, snap and subfind ID of
% all objects. 

outStruct=struct('count',[],'Simulation',[],'Snapshot',[],'SubfindID',[],...
    'ScoreRandom',[],'ScorePreferred',[]);

outStruct.count=height(objectTableComp);

outStruct.Simulation=str2double(objectTableComp.sim.extractAfter("TNG"));
outStruct.Snapshot=objectTableComp.snap;
outStruct.SubfindID=objectTableComp.subfind;
outStruct.ScoreRandom=objectTableComp.scoreRand;
outStruct.ScorePreferred=objectTableComp.scorePref;

%% write to catalog
catName='cosmic_jellyfish_viewingAngle_comparison';

folder='jellyfish_flags';

illustris.utils.write_catalog(outStruct,00,'fullName',catName,...
    'path','default','folder',folder,'v');
fprintf('===============\n');




% snaps=unique(objectTable.snap);
% sims=unique(objectTable.sim);
% 
% % run over simulations
% for k=1:length(sims)
%     bp=illustris.set_env(sims(k),'nomount');
%     global simDisplayName
%     fprintf('running on simulation %s \n ', sims(k));
%     
%     simMask=strcmp(objectTable.sim,sims(k));
%     
%     %     switch sims(k)
%     %         case 'TNG100'
%     %             funky=funkyTNG100;
%     %         case 'TNG50'
%     %             funky=funkyTNG50;
%     %     end
%     
%     fprintf('going over %i snapshots \n',length(snaps));
%     
%     % run over snap
%     
%     for i=1:length(snaps)
%         
%         fprintf('%i running on snapshot %s \n ',i, num2str(snaps(i)));
%         
%         % find the objects from the right snapshot
%         snapMask=objectTable.snap==snaps(i);
%         
%         fullMask=snapMask & simMask;
%         
%         if ~any(fullMask)  % skip this snap if there are no objects in the sim/snap
%             continue
%         end
%         
%         
%         
%         outStruct=struct('done',[],'Score',[],'ClassificationNum',[],'ScoreTotal',[]);
%         
%             
%         
%         subs=illustris.groupcat.loadSubhalos(bp,snaps(i),{'SubhaloFlag'});
%         
%         galScore=zeros(length(subs),1)-1; % set all to -1
%         galNcls=zeros(length(subs),1)-1;
%         galScoreTot=zeros(length(subs),1)-1;
%         
%         
%         
%         tabInd=find(fullMask); % indices in objectTable from the right simulation and snpashot
%         
%         % find their subfind index in the table
%         subfindInd=objectTable.subfind(fullMask)+1; % subfind is zero based, so index is increases by 1
%         
%         %prepare the scores
%         galScore(subfindInd)=objectTable.score(tabInd);
%         galNcls(subfindInd)=objectTable.clsNum(tabInd);
%         galScoreTot(subfindInd)=objectTable.scoreTotal(tabInd);
%         
%         %         %% address funky galaxyis
%         %         funkyID=funky.(['snap' num2str(snaps(i))]);
%         %
%         %         if ~isempty(funkyID)
%         %             galScore(funkyID+1)=0;
%         %             galNcls(funkyID+1)=0;
%         %             galScoreTot(funkyID+1)=0;
%         %         end
%         
%         outStruct.Score=galScore;
%         outStruct.ClassificationNum=galNcls;
%         outStruct.ScoreTotal=galScoreTot;
%         outStruct.done=int8(galScore~=-1);
%         
%         
%         
%         %% write to catalog
%         catName='cosmic_jellyfish_flags';
%         
%         folder=['jellyfish_flags/' simDisplayName];
%         
%         illustris.utils.write_catalog(outStruct,snaps(i),'name',catName,...
%             'path','default','folder',folder,'v');
%         fprintf('===============\n');
%     end
%     
% end
% 
% 
% 
% 
% 
% 
% 
% 
