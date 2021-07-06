function res = current_function()
%CURRENT_FUNCTION return the namne of the current function as a string

dbk=dbstack(1); % get current function name
if isempty(dbk)
    res="base";
else
    res=string(dbk(1).name);
end

end

