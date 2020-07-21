
function res = chooseClusterProjection
%CHOOSECLUSTERPROJECTION Interactive dialog box for choosing a cluster and a projection 
%   to plot and analyze

d = dialog('Position',[300 300 350 250],'Name','Select Cluster and Projection');
txt1 = uicontrol('Parent',d,...
    'Style','text',...
    'Position',[20 180 210 40],...
    'String','Select a cluster');

popup = uicontrol('Parent',d,...
    'Style','popup',...
    'Position',[100 170 100 25],...
    'String',{'CL101';'CL102';'CL103';...
    'CL104';'CL105';'CL106';'CL107';...
    'CL3';'CL5';'CL6';'CL7';'CL9';...
    'CL10';'CL11';'CL14';'CL24'},...
    'Callback',@popup_callback1);

txt2 = uicontrol('Parent',d,...
    'Style','text',...
    'Position',[20 120 210 40],...
    'String','Select a projection');

popup2 = uicontrol('Parent',d,...
    'Style','popup',...
    'Position',[100 100 100 25],...
    'String',{'XY','YZ','XZ'},...
    'Callback',@popup_callback2);



btn = uicontrol('Parent',d,...
    'Position',[89 20 70 25],...
    'String','Submit',...
    'Callback','delete(gcf)');

%choice = 'CL101';

% Wait for d to close before running to completion
uiwait(d);

    function popup_callback1(popup,~)
        idx = popup.Value;
        popup_items = popup.String;
        % This code uses dot notation to get properties.
        % Dot notation runs in R2014b and later.
        % For R2014a and earlier:
        % idx = get(popup,'Value');
        % popup_items = get(popup,'String');
        res.cluster = char(popup_items(idx,:));
    end

    function popup_callback2(popup,~)
        idx = popup.Value;
        popup_items = popup.String;
        % This code uses dot notation to get properties.
        % Dot notation runs in R2014b and later.
        % For R2014a and earlier:
        % idx = get(popup,'Value');
        % popup_items = get(popup,'String');
        res.proj = char(popup_items(idx,:));
    end


end





