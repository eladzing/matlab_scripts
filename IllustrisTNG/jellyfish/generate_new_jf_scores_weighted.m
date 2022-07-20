%% generate a new score table for the CJF sample based on the following:
% Incorporating an Inspector weighting scheme, which takes into account
% classification history, but also sets higher weight to known experts
% In addition, we add the classifications of the TNG50 pilot classification
% project and the Yun et al. 2019 clssifications as expert classifications.


global illUnits
global DEFAULT_MATFILE_DIR
global DEFAULT_PRINTOUT_DIR


%% read data
if readFlag
    
    % read in original CJF classification
    load([DEFAULT_MATFILE_DIR '/cosmic_jellyfish_objectTable.mat'])
    objectCJF=objectTable;
    
    objectCJF.score(objectCJF.score>20)=20;
    objectCJF.score=objectCJF.score./20;  % normalize to 1.
    clear objectTable
    
    % read in classification data for the TNG50 Pilot
    load([DEFAULT_MATFILE_DIR '/jf_objectTable_TNG50.mat'])
    load([DEFAULT_MATFILE_DIR '/jellyfishScores_TNG50.mat'])
    objectTNG50=objectTable;
    clear objectTable
    
    % read in classification data for the Yun+209 project
    load([DEFAULT_MATFILE_DIR '/yun19_jellyfish_objectTable.mat'])
    object100=objectTable;
    clear objectTable
    
    % read in Inspector score
    load([DEFAULT_MATFILE_DIR '/CJF_inspector_scores.mat']);
    
    % read in classification data for CJF
    
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
    %[idList, ia,~]=unique(clsTab.subject_ids);
    
end

%% Identify objects (object index) from CJF that are found in the previous projects.
if setupFlag
    % TNG50 pilot
    mask50=(objectCJF.sim=="TNG50");
    indCJF50=zeros(height(objectTNG50),1)-1; % tells you the index in CJF of the pilot objects
    normTNG50=6;
    for j=1:height(objectTNG50)
        subID=objectTNG50.subfind(j);
        snp=objectTNG50.snap(j);
        
        ind=find(objectCJF.snap==snp & objectCJF.subfind==subID & mask50);
        if ~isempty(ind)
            indCJF50(j)=ind;
        end
    end
    indMask50=indCJF50>0;
    score50=tally6(indMask50);%./normTNG50;  % normalized score
    
    
    
    % TNG100 Yun+2019
    mask100=(objectCJF.sim=="TNG100");
    indCJF100=zeros(height(object100),1)-1; % indices in CJF object table of matched objects
    normTNG100=5;
    for j=1:height(object100)
        subID=object100.subfind(j);
        snp=object100.snap(j);
        
        ind=find(objectCJF.snap==snp & objectCJF.subfind==subID & mask100);
        if ~isempty(ind)
            indCJF100(j)=ind;
        end
    end
    
    indMask100=indCJF100>0;
    
    %scoreCJF100=min(objectCJF.score(indCJF100(indMask100))./20,1);
    score100=double(object100.score(indMask100))./normTNG100;
end
%% generate Inspector weighting


% experts
% identify experts
iexp(1)=find(iscore.user_name.contains('zinger'));
iexp(2)=find(iscore.user_name.contains('illepic'));
iexp(3)=find(iscore.user_name.contains('rohr'));
iexp(4)=find(iscore.user_name.contains('gandha'));
iexp(5)=find(iscore.user_name.contains('dnels'));
expName=iscore.user_name(iexp);
expWeight=5;

% objectScore-down voting

strikeLimit=3;

downVoteWeight=ones(size(iscoreHigh.score));
downMask=badScore.downVote>=strikeLimit & badScore.downVoteScore>0.5;
downVoteWeight(downMask)=0;

upVoteWeight=ones(size(iscoreHigh.score));
upMask=badScore.upVote>=strikeLimit & badScore.upVoteScore>0.5;
upVoteWeight(upMask)=0;
fprintf(' %i strikes \n',strikeLimit);
fprintf('down votes %i %f \n',sum(downVoteWeight<1),sum(downVoteWeight<1)/6495);
fprintf('up votes %i %f \n',sum(upVoteWeight<1),sum(upVoteWeight<1)/6495);
fprintf('both %i %f \n' ,sum(downVoteWeight<1 & upVoteWeight<1),sum(downVoteWeight<1 & upVoteWeight<1)/6495);



wt_hi=iscoreHigh.score;
wt_hi(iscoreHigh.ncls<10)=0.5;
wt_hi_down=wt_hi.*upVoteWeight;
wt_hi_down_up=wt_hi.*upVoteWeight.*downVoteWeight;
wt_down_up=upVoteWeight.*downVoteWeight;

wt_hi_exp=wt_hi;
wt_hi_down_exp=wt_hi_down;
wt_hi_down_up_exp=wt_hi_down_up;
wt_down_up_exp=wt_down_up;

wt_hi_exp(iexp)=expWeight;
wt_hi_down_exp(iexp)=expWeight;
wt_hi_down_up_exp(iexp)=expWeight;
wt_down_up_exp(iexp)=expWeight;

%% Run over all random orientatoins  objects and generate new score
prcStep=10;
prcLvl=0;
idList=objectCJF.subject_ids;
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
    newScore.baseScore(i)=sum(votes)./length(votes);
    
    
    % identify Inspectors of this object
    unames=clsTab.user_name(inds);
    
    
    %% identify if it is in a previous study
    wt0=0;
    score0=0;
    i100=find(i==indCJF100);
    if ~isempty(i100)
        wt0=normTNG100.*expWeight;
        score0=double(object100.score(i100).*expWeight);
    end
    
    i50=find(i==indCJF50);
    
    if ~isempty(i50)
        wt0=normTNG50.*expWeight;
        score0=double(max(tally6(i50)).*expWeight);
    end
    
    %% generate new scores
    % find their scores
    [~,ii]=ismember(unames,iscore.user_name);
        
    newScore.score_hi(i)=(score0+sum(wt_hi(ii).*votes))/(wt0+sum(wt_hi(ii)));
    newScore.score_hi_down(i)=(score0+sum(wt_hi_down(ii).*votes))/(wt0+sum(wt_hi_down(ii)));
    newScore.score_hi_down_up(i)=(score0+sum(wt_hi_down_up(ii).*votes))/(wt0+sum(wt_hi_down_up(ii)));
    
    newScore.score_down_up(i)=(score0+sum(wt_down_up(ii).*votes))/(wt0+sum(wt_down_up(ii)));
    
    newScore.score_hi_exp(i)=(score0+sum(wt_hi_exp(ii).*votes))/(wt0+sum(wt_hi_exp(ii)));
    newScore.score_hi_down_exp(i)=(score0+sum(wt_hi_down_exp(ii).*votes))/(wt0+sum(wt_hi_down_exp(ii)));
    newScore.score_hi_down_up_exp(i)=(score0+sum(wt_hi_down_up_exp(ii).*votes))/(wt0+sum(wt_hi_down_up_exp(ii)));
    
    newScore.score_down_up_exp(i)=(score0+sum(wt_down_up_exp(ii).*votes))/(wt0+sum(wt_down_up_exp(ii)));
    
    %% 
    newScore.weight0(i)=wt0;
    newScore.ncls0(i)=wt0./expWeight;
    newScore.weight(i)=wt0+sum(wt_hi_down_up_exp(ii));
    newScore.ncls(i)=sum(wt_hi_down_up_exp(ii)~=0);
    newScore.nexp(i)=sum(wt_hi_down_up_exp(ii)==expWeight);    
    newScore.nclsOriginal(i)=Ncls;
end


%% repeat for the preferred direction 

prcStep=10;
prcLvl=0;
idList=objectTablePref.subject_ids;
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
    %objTab.clsNum(i)=Ncls; % total number of classifications
    
    votes=clsTab.isJelly(inds);
    newScorePref.baseScore(i)=sum(votes)./length(votes);
    
    
    % identify Inspectors of this object
    unames=clsTab.user_name(inds);
    
    
    %% identify if it is in a previous study
    wt0=0;
    score0=0;

    
    %% generate new scores
    % find their scores
    [~,ii]=ismember(unames,iscore.user_name);
        
    newScorePref.score_hi(i)=(score0+sum(wt_hi(ii).*votes))/(wt0+sum(wt_hi(ii)));
    newScorePref.score_hi_down(i)=(score0+sum(wt_hi_down(ii).*votes))/(wt0+sum(wt_hi_down(ii)));
    newScorePref.score_hi_down_up(i)=(score0+sum(wt_hi_down_up(ii).*votes))/(wt0+sum(wt_hi_down_up(ii)));
    
    newScorePref.score_down_up(i)=(score0+sum(wt_down_up(ii).*votes))/(wt0+sum(wt_down_up(ii)));
    
    newScorePref.score_hi_exp(i)=(score0+sum(wt_hi_exp(ii).*votes))/(wt0+sum(wt_hi_exp(ii)));
    newScorePref.score_hi_down_exp(i)=(score0+sum(wt_hi_down_exp(ii).*votes))/(wt0+sum(wt_hi_down_exp(ii)));
    newScorePref.score_hi_down_up_exp(i)=(score0+sum(wt_hi_down_up_exp(ii).*votes))/(wt0+sum(wt_hi_down_up_exp(ii)));
    
    newScorePref.score_down_up_exp(i)=(score0+sum(wt_down_up_exp(ii).*votes))/(wt0+sum(wt_down_up_exp(ii)));
    
    %% 
%     newScorePref.weight0(i)=wt0;
%     newScorePref.ncls0(i)=wt0./expWeight;
    newScorePref.weight(i)=wt0+sum(wt_hi_down_up_exp(ii));
    newScorePref.ncls(i)=sum(wt_hi_down_up_exp(ii)~=0);
    newScorePref.nexp(i)=sum(wt_hi_down_up_exp(ii)==expWeight);    
    newScorePref.nclsOriginal(i)=Ncls;
end






%% save new score table

