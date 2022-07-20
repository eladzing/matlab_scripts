if setupFlag
    %% read in table and format it to our liking.
    clsTab=jellyfish.readClassification.import_classification_file_CJF('191121');
    %% separate the Phase 1 classification table
    
    startDate=datetime('01/06/21 00:00','InputFormat',"dd/MM/yy HH:mm");
    endDate=datetime('14/06/21 23:59','InputFormat',"dd/MM/yy HH:mm");
    
    clsTabP1=jellyfish.readClassification.format_classification_table_CJF(clsTab,'start',startDate,'end',endDate);
    %head(clsTab,5)
    %clsTab=jellyfish.readClassification.simulationMask(clsTab,'TNG50');
    % Classification info based on data downloaded on 14/7/21
    % review and training set objecs classifications are removed.
    % The raw classification data is reformatted into a new table with columns contiaing the snapshot, id  Id for the object as well as the label (jelly or not):
    %% Separaete the phase 2 classification table
    
    startDate=datetime('22/08/21 00:00','InputFormat',"dd/MM/yy HH:mm");
    endDate=datetime('19/11/21 23:59','InputFormat',"dd/MM/yy HH:mm");
    
    clsTabP2=jellyfish.readClassification.format_classification_table_CJF(clsTab,'start',startDate,'end',endDate);
    
    %% concatenate the two phases
    
    clsTab12=vertcat(clsTabP1,clsTabP2);
    GC=groupcounts(clsTab12,'user_name');
    
    clear clsTabP1 clsTabP2 clsTab
    %% load object table
    
    global DEFAULT_MATFILE_DIR
    load([DEFAULT_MATFILE_DIR '/cosmic_jellyfish_objectTable.mat'],'objectTable','objectTablePref');
    
    sidList=cat(1,objectTable.subject_ids(:),objectTablePref.subject_ids(:));
    clear objectTable objectTablePref
    
end
%% Generate an Inspector score
% We will generate an inspector score - based on how often the inspector votes
% differently than the majority. The higher/lower the overall score, the more
% the discrepancy will count. We'll do this for high and low scores separately
% and also for the combined



% goodness score for Inspectors
iscore=GC;
iscore.wt=zeros(size(iscore.user_name));
iscore.score=zeros(size(iscore.user_name));
iscore.ncls=zeros(size(iscore.user_name));

iscoreHigh=iscore; % inspector score based on voting record in high-score objects
iscoreLow=iscore;  % inspector score based on voting record in low-score objects

% identifying down-voters and up-voters
badScore=GC;
badScore.upVote=zeros(size(badScore.user_name));
badScore.upVoteScore=zeros(size(badScore.user_name));
badScore.upVoteNum=zeros(size(badScore.user_name));
badScore.upVoteWt=zeros(size(badScore.user_name));

badScore.downVote=zeros(size(badScore.user_name));
badScore.downVoteScore=zeros(size(badScore.user_name));
badScore.downVoteNum=zeros(size(badScore.user_name));
badScore.downVoteWt=zeros(size(badScore.user_name));

vFrac=0.75:0.05:1;
%upVFrac=0:0.05:0.25;       1-downVFrac;


prcStep=10;
prcLvl=0;
fprintf('Running over %i objects \n \n',length(sidList))
for i=1:length(sidList)
    
    prc=round(i./length(sidList).*100);
    if prc>=prcLvl
        fprintf('completed %i %% of objects \n',prc);
        prcLvl=prcLvl+prcStep;
    end
    
    
    clsInd=find(clsTab12.subject_ids==sidList(i)); % indices of classifications
    
    unames=clsTab12.user_name(clsInd);  % names of classifiers
    votes=clsTab12.isJelly(clsInd);   % votes of Ispectors
    objectScore=sum(votes)./length(votes);
    [~,ii]=ismember(unames,GC.user_name);  % indices in the table
    
    wt0=2.0.*objectScore-1;
    wt=abs(wt0);
    
    vv=votes;
    
    
    if sign(wt0)>0
        iscoreHigh.wt(ii)=iscoreHigh.wt(ii)+wt;
        iscoreHigh.score(ii)=iscoreHigh.score(ii)+wt.*vv;
        iscoreHigh.ncls(ii)=iscoreHigh.ncls(ii)+1;
    elseif sign(wt0)<0
        vv=1-vv;
        
        iscoreLow.wt(ii)=iscoreLow.wt(ii)+wt;
        iscoreLow.score(ii)=iscoreLow.score(ii)+wt.*vv;
        iscoreLow.ncls(ii)=iscoreLow.ncls(ii)+1;      
    end
    
    
    iscore.wt(ii)=iscore.wt(ii)+wt;
    iscore.score(ii)=iscore.score(ii)+wt.*vv;
    iscore.ncls(ii)=iscore.ncls(ii)+1;
    
    
    %% upvotes and down votes
    
    % find wieght according to score
    
    wtHi=discretize(objectScore,vFrac);
    wtLo=discretize(1-objectScore,vFrac);
    
    if ~isnan(wtLo) % score is in the low range
        badScore.upVote(ii)=badScore.upVote(ii)+votes;
        badScore.upVoteScore(ii)=badScore.upVoteScore(ii)+wtLo.*votes;
        badScore.upVoteNum(ii)=badScore.upVoteNum(ii)+1;
        badScore.upVoteWt(ii)=badScore.upVoteWt(ii)+wtLo;
    elseif ~isnan(wtHi) % score is in the low range
        badScore.downVote(ii)=badScore.downVote(ii)+~votes;
        badScore.downVoteScore(ii)=badScore.downVoteScore(ii)+wtHi.*(1-votes);
        badScore.downVoteNum(ii)=badScore.downVoteNum(ii)+1;
        badScore.downVoteWt(ii)=badScore.downVoteWt(ii)+wtHi;
%     else  % score is in the middle
%         continue
    end
    
end

iscore.score=iscore.score./iscore.wt;

iscore.score(iscore.wt==0)=1;



iscoreHigh.score=iscoreHigh.score./iscoreHigh.wt;
iscoreHigh.score(iscoreHigh.ncls==0)=1;

iscoreLow.score=iscoreLow.score./iscoreLow.wt;
iscoreLow.score(iscoreLow.ncls==0)=1;

upmask=badScore.upVoteNum>0;
badScore.upVoteScore(upmask)=badScore.upVoteScore(upmask)./badScore.upVoteWt(upmask);
downmask=badScore.downVoteNum>0;
badScore.downVoteScore(downmask)=badScore.downVoteScore(downmask)./badScore.downVoteWt(downmask);





%% save to mat file
save([DEFAULT_MATFILE_DIR '/CJF_inspector_scores.mat'],...
    'iscore','iscoreHigh','iscoreLow','badScore');