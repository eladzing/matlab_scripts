function res = fixstring( strin )
%FIXSTRING fixes underscores in strings for use in latex.

if ~ischar(strin)
    error('fixstring: argument must be a string')
end

ki=strfind(strin,'_');

if ~isempty(ki)
    ind=1;
    st=cell(length(ki)+1,1);
    for i=1:length(ki)
        st{i}=strcat(strin(ind:ki(i)-1),'\_');
        ind=ki(i)+1;
    end
    st{end}=strin(ind:end);
    
    res=st{1};
    for i=2:length(st)
        res=strcat(res,st{i});
    end
else
    res=strin;
end



end

