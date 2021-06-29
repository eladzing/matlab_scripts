function hand=titlemine(string,varargin)
%% my own version of title set to latex and a default fontsize of 14
% @param string is the string in the label 
% @param varargin is an optional input for size



%Don't want more than 2 arugments altogether
numvarargs=length(varargin);
if numvarargs>1
    error('titlemine: too many inputs','Only two are allowed');
end

%set defualt  fontsize
defvals={14};

%assign the optional values
defvals(1:numvarargs)=varargin;

% transfer to easy to use varaibles
size=defvals{:};

hand=title(string,'Fontsize',size,'Interpreter','latex');




