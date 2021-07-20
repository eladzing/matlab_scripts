function res = build_object_table(clsTab)
%BUILD_OBJECT_TABLE Given a table of classifications, constrct a table of
%unique objects with identifying information and scores 
%   Detailed explanation goes here

% find unique subjects
[idList, ia,~]=unique(clsTab.subject_ids);

% create object table
res=table(idList,clsTab.simName(ia),clsTab.snap(ia),clsTab.subfind(ia),zeros(size(idList)),zeros(size(idList)),...
    'variableNames',{'subject_ids','sim','snap','subfind','Ncls','Ycls'});

%loop over subject ids and count classifications
for i=1:length(idList)
    inds=find(clsTab.subject_ids==idList(i));
    
    res.Ncls(i)=length(inds); % total number of classifications
    
    res.Ycls(i)=sum(clsTab.isJelly(inds)); % number of 'yes' answers
    
    
end

end

