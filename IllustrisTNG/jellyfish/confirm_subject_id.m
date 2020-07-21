%% confirm that subject ID is one-to-one with subfind/snap   - CONFIRMED! 

% inds=find(clsTab.workflow_id==11990);
% 
% 
% for i=1:length(inds)
%     subID(i)=clsTab.subject_ids(inds(i));
%     subInfo(i)=get_subject_info(clsTab.subject_data(inds(i)));
%     
% end
%%

uID=unique(subID);

conf=true(length(uID),5);
for i=1:length(uID)
    
    ind2=find(subID==uID(i));
    
    conf2=true(length(ind2),5);
    for j=1:length(ind2)
        for k=1:5
            switch k
                case 1
                    conf2(j,k)= uID(i)==subInfo(ind2(j)).subjectID;
                    
                case 2
                    conf2(j,k)= strcmp(subInfo(ind2(1)).fullSimName,subInfo(ind2(j)).fullSimName);
                    
                case 3
                    conf2(j,k)= subInfo(ind2(1)).snap==subInfo(ind2(j)).snap;
                case 4
                    conf2(j,k)= subInfo(ind2(1)).hostID==subInfo(ind2(j)).hostID;
                case 5
                    conf2(j,k)= subInfo(ind2(1)).subfindID==subInfo(ind2(j)).subfindID;
            end
        end
    end
    
    for k=1:5
        conf(i,k)=all(conf2(:,k));
    end
    
end
    
    
    