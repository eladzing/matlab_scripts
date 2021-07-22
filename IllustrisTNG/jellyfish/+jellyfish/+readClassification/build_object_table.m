function res = build_object_table(clsTab)
%BUILD_OBJECT_TABLE Given a table of classifications, construct a table of
%unique objects with identifying information and scores
%   Detailed explanation goes here

%
clsBase=20;
Nsample=200;

% find unique subjects
[idList, ia,~]=unique(clsTab.subject_ids);

% create object table
res=table(idList,clsTab.simName(ia),clsTab.snap(ia),clsTab.subfind(ia),zeros(size(idList)),zeros(size(idList)),...
    'variableNames',{'subject_ids','sim','snap','subfind','clsNum','score','scoreTotal'});

%loop over subject ids and count classifications
for i=1:length(idList)
    inds=find(clsTab.subject_ids==idList(i));
    
    Ncls=length(inds);
    res.clsNum(i)=Ncls; % total number of classifications
    
    isJel=clsTab.isJelly(inds);
    res.scoreTotal(i)=sum(isJel); % number of 'yes' answers
    
    % deal with objects with more than the base number of classifications
    if Ncls<=clsBase
        res.score(i)=sum(isJel);
    else
        for j=1:Nsample
            [~,ix]=sort(rand(1,Ncls));
            scores(j)=sum(isJel(ix(1:clsBase)));
        end
        res.score(i)=median(scores);
    end
    
    
end

end

