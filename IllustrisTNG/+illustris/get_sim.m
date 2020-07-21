function simulation = get_sim(simName)
%ENV_API Set basic environment for Illustris data acquistion through API
%   Detailed explanation goes here

% Get info about simulations: name, # of snapshots and url
%simInfoStruct = get_url(baseUrl);

% Extract simulation names
names={};
for i=1:numel(simInfoStruct.('simulations'))
    names{i} = simInfoStruct.('simulations'){i}.('name');
end

if ~exist('simName','var')
    
    % get data interactively
    %simName=uicontrol('Style','popup','string',names,'position',[])
    
    % dh = dialog('Position',[800 800 300 100],'Name','Choose Simulation');
    simName=choosedialog(names);%   uicontrol('parent',dh,'Style','popup','string',names,'Position',[75 70 100 25],'Callback',@popup_callback);
end

[~,simInd] = ismember(simName,names);

simulation = get_url(simInfoStruct.('simulations'){simInd}.('url') );
end

function choice = choosedialog(names)

%choice='';%names{1};
dh = dialog('Position',[800 500 300 200],'Name','Choose Simulation');
ch='';
uicontrol('Parent',dh,...
    'Style','text',...
    'units','normalized',...
    'Position',[0.2 0.5 0.4 0.3],...
    'String','Select a simulation');
popup=uicontrol('parent',dh,'Style','popup','string',names,...
    'units','normalized','Position',[0.2 0.2 0.5 0.3],...
    'Callback',@popup_callback);

btn = uicontrol('Parent',dh,...
    'Position',[89 20 70 25],...
    'String','Close',...
    'Callback','delete(gcf)');

%choice=popup.String{popup.Value};
 uiwait(dh);
choice=ch;


    function popup_callback(pop,~)
        %      idx = popup.Value;
        %           popup_items = popup.String;
        ch=pop.String{pop.Value};
        %delete(gcf)
        
    end

end