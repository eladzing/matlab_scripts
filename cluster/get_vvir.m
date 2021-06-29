function result=get_vvir(varargin)
%% get vvir/v200/v500 from list.
%optional argument can be empty or 200 or 500

if isempty(varargin)
    res=get_virials();
else
    d=varargin{:};
    res=get_virials(d);
end
    
result=res(3);

end

    
    
    
    