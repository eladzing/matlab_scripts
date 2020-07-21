function result=get_mvir(varargin)
%% get mvir/m200/m500 from list.
%optional argument can be empty or 200 or 500

if isempty(varargin)
    res=get_virials();
else
    d=varargin{:};
    res=get_virials(d);
end
    
result=res(2);

end

    
    
    
    