function matched=create_matched_sample(groupID,poolID,param)
%% description 


matched=zeros(size(groupID));

tempList=poolID;
    

[~, newInd]=sort(rand(size(groupID)));
    
    groupID2=groupID(newInd);
    
    
    for i=1:length(groupID2)
        
        dist=abs(galMass(groupID2(i))-galMass(tempList));
        [~,choice]=min(dist);
        
        cand=tempList(choice);
        matched(i)=cand;
        
        tempList=tempList(tempList~=cand);
        
    end
    
    
end
%%
%figure;loglog(galMass(groupID2),galMass(matched),'.')

%figure;semilogx(galMass(groupID2),galMass(matched)./galMass(groupID2)-1,'.')