function barTitle(barHandle,titleString,varargin)
%%  function to set the title over a colorbar
%
% Syntax: 
% barTitle(barHandle,titleString) - sets the title of the bar with handle
% barHandle to be titleString. The defualt interpreter os latex and the
% defualt Font is 14. 
%
% optional arguments:
% interpreter 
% fontsize

%% defualts 
fontSize=14;
interpret='latex';
col='k';
i=1;
while i<=length(varargin)
    switch lower(varargin{i})
        case 'Interpreter'
            i=i+1;
            interpret=varargin{i};
        case {'fontsize','font','fontSize','FontSize','Fontsize'}
            i=i+1;    
            fontSize=varargin{i};
        case{'color','col','fontcolor'}
            i=i+1;
            col=varargin{i};
        otherwise
            error('barTitle: Illegal option %s',varargin{i})
    end
    i=i+1;
end


set(get(barHandle,'Title'),'String',titleString,'Fontsize',fontSize,'Interpreter',interpret,'color',col)