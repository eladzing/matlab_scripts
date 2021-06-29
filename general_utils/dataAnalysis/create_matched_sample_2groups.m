function res=create_matched_sample_2groups(groupID,poolID,param,method)
%% CREATE_MATCHED_SAMPLE - create a matched sample to given group
% given a group of objects, find a matched group, in some parameter space
% for a given object, the match pair is found by finding the closest object
% in the 'pool' in terms of the given parameter.
% Once a matched pair is made, the matched object is removed from the
% 'pool' to avoid doubl matching.
% The initial group is scrambled to avoid trends.
% inputs arguments: groupID and poolID are lists of indices pointing to
% values in the param array


matched=zeros(size(groupID));
tempList=poolID;

%% mix up the original sample
aa=rand(size(groupID));
[~, newInd]=sort(aa);
groupID2=groupID(newInd);

oldInd(newInd) =1:length(aa);

if ~exist('method','var')
    method='minimal';
end


switch(lower(method))
    
    case{'minimal','min','minmatch'}
        
        for i=1:length(groupID2)
            
            % find the closest match
            dist=abs(param(groupID2(i))-param(tempList));
            [~,choice]=min(dist);
            
            
            matched(i)=tempList(choice);
            
            % remove the match from the pool
            tempList=tempList(tempList~=matched(i));
            
        end
        
    case{'range','rangematch'}
        
        massFac=0.05;
        step=0.05;
        add=0;
        
        i=1;
        while i<=length(groupID2)
            
            range=param(groupID2(i)).*(1+[-1 1].*(massFac+add));
            
            mm=param(tempList)>=range(1) & param(tempList)<=range(2);
            
            if sum(mm)>1
                
                cand=tempList(mm);
                
                ran=rand(size(cand));
                
                [~,choice]=max(ran);
                
                matched(i)=cand(choice);
                
                tempList=tempList(tempList~=cand(choice));
                
                i=i+1;
            else
                add=add+step;
            end
            
            if massFac>0.5
                
                error('CREATE_MATCHED_SAMPLE - range exceeds 50%')
            end
        end
        
    otherwise
        error('CREATE_MATCHED_SAMPLE - Illegal method: %s',method)
end




res=matched(oldInd);

end
%%
%figure;loglog(param(groupID2),param(matched),'.')

%figure;semilogx(param(groupID2),param(matched)./param(groupID2)-1,'.')