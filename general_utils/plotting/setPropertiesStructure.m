function setPropertiesStructure(handle,propStruct)
%SETPROPERTIESSTRUCTURE helper function to convert a structure which
%contans properties of an  object to a SET command for these properties for
%a given handle

fields=fieldnames(propStruct);
for i=1:length(fields)
    
    set(handle,fields{i},propStruct.(fields{i}))
    
end

end

