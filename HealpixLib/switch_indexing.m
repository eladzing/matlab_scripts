function indNew=switch_indexing(rindex)
%% using the ring indexing, created an index array to arrange the cells in the Ring scheme
ind=1:size(rindex,1);

for i=1:max(rindex(:,1))
    ff1=find(rindex(:,1)==i);
    [~,ff2]=sort(rindex(ff1,2));
    
    aa=ind(ff1);
    bb=aa(ff2);
    
    if i==1
        indNew=bb;
    else
        indNew=cat(2,indNew,bb);
    end
end

end
    
