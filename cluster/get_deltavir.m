function result=get_deltavir(varargin)
%% get the delta_vir used to calculate the virials  from list.
%optional argument can be empty or 200 or 500

if isempty(varargin)
    res=get_virials();
else
    d=varargin{:};
    res=get_virials(d);
end
    
result=res(5);

end

    
    
    
    