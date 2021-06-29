function result = Tgal(aexpn)
% Utility function to load the gas temperature cube
%
% @param aexpn The required expansion parameters
% @returns  A 300^3 matrix.
%

narginchk(1,1);

%aexp=cell(size(aexpn));


global GAL_PATH;
if ischar(aexpn)
    aexp{1}=aexpn;
else
    aexp=aexpn;
end

for i=1:length(aexp)
    fnameBase=[GAL_PATH '/a' aexp{i} '*_temperature.bin'];
    
    fname=ls(fnameBase);
    
    if isempty(fname)
        error('Tgal - Nonexistant expansion parameter: %s',aexp{i})
    end
    
    fnameSize=size(fname);
    
    if (fnameSize>1)
        error('Tgal - More than one file fits the bill')
    end
    
    filename=[GAL_PATH '/' fname];
    result(i)= read_galCube(filename);
    a_exp(i)=str2double(aexp{i});
end

for i=1:length(a_exp)
    result(i).aexp=a_exp(i);
    result(i).zred=1/a_exp(i)-1;
end


end
