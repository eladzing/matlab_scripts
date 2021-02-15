for k=1:3
    listMask=false(size(satStruct(k).satID));
    
    for j=1:length(satStruct(k).satID)
        
        indlist=find(satStruct(k).satID==satStruct(k).satID(j));
        if length(indlist)==1
            listMask(indlist)=true;
        else
            ps=satStruct(k).pscore(indlist);
            [~,ix]=max(ps);
            listMask(indlist(ix))=true;
            
        end
        
    end
    
    satStruct(k).listMask=listMask;
end