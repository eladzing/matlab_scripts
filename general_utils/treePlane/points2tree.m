function cellVal = points2tree(val,treeStruct,method)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~exist('method','var')
    method='avg';
end

cellVal=zeros(size(treeStruct.treeLevel));

for i=1:length(treeStruct.treeLevel)
    
    maskP=treeStruct.treeBLC(1,i)==treeStruct.pointsBLC(1,:) & ...
        treeStruct.treeBLC(2,i)==treeStruct.pointsBLC(2,:);
    
    
    vv=val(maskP);
    mask=~isnan(vv) & ~isinf(vv);
    
    
    switch(lower(method))
        case{'avg','mean'}
            
            cellVal(i)=sum(vv(mask))./sum(mask);
            
        case{'sum','total'}
            
            cellVal(i)=sum(vv(mask));
            
        case{'median','med'}
            cellVal(i)=median(vv(mask));
        otherwise
            error('points2tree - Illegal method: %s',method)
    end
end



end

