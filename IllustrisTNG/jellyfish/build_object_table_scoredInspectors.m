% function res = build_object_table(clsTab)
% %BUILD_OBJECT_TABLE Given a table of classifications, construct a table of
% %unique objects with identifying information and scores
% %   Detailed explanation goes here

%
%clsBase=20;
%Nsample=200;
if setupFlag
    %% create classification table
    clsTab0=jellyfish.readClassification.import_classification_file_CJF('191121');
    %% separate the Phase 1 classification table
    startDate=datetime('01/06/21 00:00','InputFormat',"dd/MM/yy HH:mm");
    endDate=datetime('14/06/21 23:59','InputFormat',"dd/MM/yy HH:mm");
    
    clsTabP1=jellyfish.readClassification.format_classification_table_CJF(clsTab0,'start',startDate,'end',endDate);
    
    % Classification info based on data downloaded on 14/7/21
    % review and training set objecs classifications are removed.
    % The raw classification data is reformatted into a new table with columns contiaing the snapshot, id  Id for the object as well as the label (jelly or not):
    
    %% Separaete the phase 2 classification table
    
    startDate=datetime('22/08/21 00:00','InputFormat',"dd/MM/yy HH:mm");
    endDate=datetime('19/11/21 23:59','InputFormat',"dd/MM/yy HH:mm");
    
    clsTabP2=jellyfish.readClassification.format_classification_table_CJF(clsTab0,'start',startDate,'end',endDate);
    %% concatenate the two phases
    
    clsTab=vertcat(clsTabP1,clsTabP2);
    
    clear clsTabP1 clsTabP2 clsTab0
    
    %% load Inspector scores
    
    
    global DEFAULT_MATFILE_DIR
    load([DEFAULT_MATFILE_DIR '/CJF_inspector_scores.mat']);
    
    
end
% find unique subjects
[idList, ia,~]=unique(clsTab.subject_ids);

% create object table
objTab=table(idList,clsTab.simName(ia),clsTab.snap(ia),clsTab.subfind(ia),clsTab.tag(ia),...
    zeros(size(idList)),...
    'variableNames',{'subject_ids','sim','snap','subfind','tag','clsNum'});

scores.scoreBase=zeros(size(idList));
scores.scoreHigh=zeros(size(idList));
scores.scoreLow=zeros(size(idList));
scores.scoreOverall=zeros(size(idList));
scores.scoreCombined=zeros(size(idList));
scores.noUpVote=zeros(size(idList));
scores.noDownVote=zeros(size(idList));
scores.noBad=zeros(size(idList));
scores.experts5=zeros(size(idList));
scores.scoreExperts5=zeros(size(idList));
scores.scoreCombinedExperts5=zeros(size(idList));
scores.scoreAll=zeros(size(idList));

upVoters=badScore.user_name(badScore.upVoteScore>0.5);
downVoters=badScore.user_name(badScore.downVoteScore>0.5);

% identify experts
iexp(1)=find(iscore.user_name.contains('zinger'));
iexp(2)=find(iscore.user_name.contains('illepic'));
iexp(3)=find(iscore.user_name.contains('rohr'));
iexp(4)=find(iscore.user_name.contains('gandha'));
iexp(5)=find(iscore.user_name.contains('dnels'));

expName=iscore.user_name(iexp);


%loop over subject ids and count classifications

prcStep=10;
prcLvl=0;
fprintf('Running over %i objects \n \n',length(idList))
for i=1:length(idList)
    
    % advancement ticker
    prc=round(i./length(idList).*100);
    if prc>=prcLvl
        fprintf('completed %i %% of objects \n',prc);
        prcLvl=prcLvl+prcStep;
    end
    
    % find all classifications     
    inds=find(clsTab.subject_ids==idList(i));
    
  
    
    Ncls=length(inds); 
    objTab.clsNum(i)=Ncls; % total number of classifications
    
    votes=clsTab.isJelly(inds); 
    
    scores.scoreBase(i)=sum(votes)./length(votes); % fraction of 'yes' answers
    
       
    
    % identify Inspectors of this object
    unames=clsTab.user_name(inds);
    
    % find their scores
    [~,ii]=ismember(unames,iscore.user_name);
    
    sc=iscore.score(ii);
    scores.scoreOverall(i)=sum(votes.*sc)./sum(sc);
    
    sc=iscoreHigh.score(ii);
    scores.scoreHigh(i)=sum(votes.*sc)./sum(sc);
    
    sc=iscoreLow.score(ii);
    scores.scoreLow(i)=sum(votes.*sc)./sum(sc);
    
    sc=0.5.*(iscoreLow.score(ii) + iscoreHigh.score(ii));
    scores.scoreCombined(i)=sum(votes.*sc)./sum(sc);
    
    % bad scores - up/down voters are ignored
    [~,iup]=ismember(unames,upVoters);
    [~,idn]=ismember(unames,downVoters);
    [~,iex]=ismember(unames,expName);
    
    sc=ones(size(votes));
    sc(iup~=0)=0;
    scores.noUpVote(i)=sum(votes.*sc)./sum(sc);
    
    sc=ones(size(votes));
    sc(idn~=0)=0;
    scores.noDownVote(i)=sum(votes.*sc)./sum(sc);
    
    sc=ones(size(votes));
    sc(idn~=0)=0;
    sc(iup~=0)=0;
    scores.noBad(i)=sum(votes.*sc)./sum(sc);
    
    sc=ones(size(votes));
    sc(iex~=0)=5;
    scores.experts5(i)=sum(votes.*sc)./sum(sc);
    
    sc=iscore.score(ii);
    sc(iex~=0)=5;
    scores.scoreExperts5(i)=sum(votes.*sc)./sum(sc);
    
    sc=0.5.*(iscoreLow.score(ii) + iscoreHigh.score(ii));
   sc(iex~=0)=5;
    scores.scoreCombinedExperts5(i)=sum(votes.*sc)./sum(sc);
    
    sc=0.5.*(iscoreLow.score(ii) + iscoreHigh.score(ii));
    sc(iex~=0)=5;
    sc(idn~=0)=0;
    sc(iup~=0)=0;
    scores.scoreAll(i)=sum(votes.*sc)./sum(sc);
    
    
end





