function res = mask_structure(inStr,mask)
%MASK_STRUCTURE enforce a mask on all arrays of a structure
%   recursively travel through a strucute and enforce a mask on every
%   dimension that has the same size as the mask
current_function='MASK_STRUCT';  % update later with real function

if ~isstruct(inStr)
    error([current_function ' - input argument not a structure'])
end

fnames=fieldnames(inStr);

len=length(mask);  % in the future we can expand this to multi-dimnsional mask
szMask=size(mask);
for i=1:length(fnames)
    
    
    if isnumeric(inStr.(fnames{i})) || islogical(inStr.(fnames{i}))
        % mask numeric or logical fields
        A=inStr.(fnames{i}); 
        sz=size(A);
        
        if sz==szMask  % mask is same shape as field then its easy
            inStr.(fnames{i})=A(mask);
        elseif any(sz==len) 
            % enforce the mask only on the relavent dimension (equal to
            % length of mass)
            S.subs=repmat({':'},1,ndims(A));
            S.type='()';
            S.subs{size(A)==len}=find(~mask);
            B=subsasgn(A,S,[]);
            inStr.(fnames{i})=B;
            
        else % if the field has no dimension equal to mask, leave it be.
            continue;
        end
    elseif isstruct(inStr.(fnames{i})) % if field is a struct - recursively change it too
        inStr.(fnames{i})=mask_structure(inStr.(fnames{i}),mask);
        
    end
end

res=inStr;
end

