function res = single_object_catalog(id,dataStr)
%SINGLE_OBJECT_CATALOG generate a structure containing the catalog values for a single fof/subhalo 
%   Detailed explanation goes here

indx=id+1;

fldNames=fieldnames(dataStr); 
for i=1:length(fldNames)
    if strcmp(fldNames{i},'count'); continue; end
    res.(fldNames{i})=dataStr.(fldNames{i})(:,indx);
end


end

