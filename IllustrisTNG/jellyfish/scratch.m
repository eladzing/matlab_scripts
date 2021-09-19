for i=1:3%length(listInd) 
    ind=listInd(i);
    snp=objectCJF.snap(ind);
    sim=objectCJF.sim(ind);
    score=objectCJF.score(ind);
    subID=objectCJF.subfind(ind);
    
    % check sim 
    if sim~="TNG100"
        error('wrong sim')
    end
    
    
    % check to see that the score is indeed 5 
    Yind=find(object100.snap==snp & object100.subfind==subID);
    if ~isempty(ind)
        Yscore=object100.score(Yind);
        if Yscore~=5
            error('wrong score')
        end
    else
        error('empty index')
    end
    
    % write to file 
    fmt='%2i  %7i  %2i  %1i \n';
    fprintf(fid,fmt,snp,subID,score,Yscore);
    
end