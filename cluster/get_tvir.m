function result=get_tvir(varargin)
%% get tvir/t200/t500 from list.
%optional argument can be empty or 200 or 500

if isempty(varargin)
    res=get_virials();
else
    d=varargin{:};
    res=get_virials(d);
end
    
result=res(4);

end

    
    
    
    