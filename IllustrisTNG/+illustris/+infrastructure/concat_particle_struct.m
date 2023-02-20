function res = concat_particle_struct(struct1,struct2)
%CONCAT_PARTICLE_STRUCTURE stitch together 2 particle structres of the same type
% if the structures contain different fields only the common fields will
% appear in the resultant structure; 

missingFlag=false;

fldName1=fieldnames(struct1);
fldName2=fieldnames(struct2);
missingFlag1=true(size(fldName1));
missingFlag2=true(size(fldName2));

fldList={};
% check if structures contain the same fields 
for i=1:length(fldName1)
    for j=i:length(fldName2)
        if strcmp(fldName1{i},fldName2{j})
            fldList{end+1}=fldName1{i};
            missingFlag1(i)=false;
            missingFlag2(j)=false;
            break
        end
    end
end

if any(missingFlag1) || any(missingFlag2)
    warning('%s - structures contain different fields. Only common fields will be considered',...
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

