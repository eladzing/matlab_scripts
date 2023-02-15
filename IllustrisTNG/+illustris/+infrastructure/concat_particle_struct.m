function res = concat_particle_struct(struct1,struct2)
%CONCAT_PARTICLE_STRUCTURE stitch together 2 particle structres of the same type
% if the structures contain different fields only the common fields will
% appear in the resultant structure; 

missingFlag=false;

fldName1=fieldnames(struct1);
fldName2=fieldnames(struct2);

fldList={};
% check if structures contain the same fields 
for i=1:length(fldName1)
    for j=1:length(fldName2)
        if strcmp(fldName1{i},fldName2{j})
            fldList{end+1}=fldName1{i};
            continue
        end
        missingFlag=true;
    end
end

if missingFlag
    warning('%s - structures contain different fields. Only common fields will be ocnsidered',...
        current_function().upper);
end

for i=1:length(fldList)
    if strcmp(fldList{i},'count')
        res.count=struct1.count+struct2.count;
    else
        res.(fldList{i})=cat(2,struct1.(fldList{i}),struct2.(fldList{i}));
    end
end


end

