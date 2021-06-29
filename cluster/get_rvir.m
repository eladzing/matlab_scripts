function result=get_rvir(varargin)
%% get rvir/r200/r500 from list.
%optional argument can be empty or 200 or 500

if isempty(varargin)
    res=get_virials();
else
    d=varargin{:};
    res=get_virials(d);
end
    
result=res(1);

end

    
    
    
    