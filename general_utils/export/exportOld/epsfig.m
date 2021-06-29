function epsfig(filename,varargin)

% Create a nice eps fig 
%  Don't want more than 3 printing arugments altogether
numvarargs=length(varargin);
if numvarargs>1
    error('epsfig: too many printing arguments')
end
    
%set defualt printing flag, printing tag, and output dir
defvals={gcf ; 'printout'};

%assign the optional values
defvals(1:numvarargs)=varargin;

% transfer to easy to use varaibles
[hfig printdir]=defvals{:};


% set sizes 
%set(hfig, 'PaperUnits', 'centimeters');
%set(hfig, 'PaperSize', [13 12]);
set(hfig, 'PaperPositionMode', 'auto');
%set(hfig, 'PaperPosition', [0 0 12 12]);

% remove whitspace 
%set(gca, 'Position', get(gca, 'OuterPosition') - ...
%    get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);

set(hfig, 'renderer', 'painters');
print(hfig, '-depsc2', sprintf('%s/%s.eps',printdir,filename));




