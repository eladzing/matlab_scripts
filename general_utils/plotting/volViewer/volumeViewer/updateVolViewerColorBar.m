% updateVolViewerColorBar - updates VolumeViewre colorbar display
%
% Usage:
% ------
% [clim, handles] = updateVolViewerColorBar(handles, colorbarStatus)
% Inputs are:
%   - handles - the object handles structure for the VolumeViewer GUI.
%   - colorbarStatus - the status string for the colorbar including the
%     following options:
%       - 'on' - turn the colorbar display on.
%       - 'off' - turn the folorbar display off.
%       - 'update' - update the colorbar display with the current settings.
% Outputs are:
%    - clim - the color limits of the colorbar.
%    - hangles - the VolumeViewre GUI object handles structure.
%
% By Ran Klein 

function [clim, handles] = updateVolViewerColorBar(handles, colorbarStatus)

if nargin<2
	colorbarStatus = 'update';
end

cmap = get(handles.volViewerFigure,'colormap');
clim = get(handles.DisplayAxis,'Clim');

if ColorbarNotInitialized(handles)
	handles = initializeColorbar(handles, clim, cmap);
else
	steps = size(cmap,1);
	set(findobj(get(handles.Colorbar,'children'),'flat','type','image'),...
		'cdata',reshape(cmap,1,steps,3),'xdata',linspace(clim(1),clim(2),steps)');
	set(handles.Colorbar,'XLim',clim,'XTickMode','Auto')
end


switch lower(colorbarStatus)
	case 'off'
		h = [handles.Colorbar; get(handles.Colorbar,'Children')];
		set(h,'Visible','off');
	case {'on','update'}
		if strcmpi(colorbarStatus,'on')
			h = [handles.Colorbar; get(handles.Colorbar,'Children')];
			set(h,'Visible','on');
		end
		% update the width of all the arrow patches so they do not become
		% distorted as the colorbar scale changes.
		if strcmp(get(handles.Colorbar,'Visible'),'on')
			h = get(handles.Colorbar,'children');
			h = findobj(h,'flat','type','patch');
			for i=1:length(h)
				pos = get(h(i),'vertices');
				pts = pos(:,1)-pos(1,1);
				pts = pts/(max(pts)-min(pts));
				scale=abs(clim(1)-clim(2))/15/1;
				pos(:,1)=pts*scale+pos(1,1);
				set(h(i),'vertices',pos);
			end
			updateContour(handles)
		end
	otherwise
		msgID = sprintf('unrecognizedOption:%s:volViewer',mfilename);
		error(msgID,['The status ' status ' is unrecognized.']);
end





function answer = ColorbarNotInitialized(handles)
answer = isempty(findobj(get(handles.Colorbar,'children'),'flat','type','image'));


function handles = initializeColorbar(handles, clim, cmap)
steps = size(cmap,1);
vol = getappdata(handles.DisplayAxis,'SummedCroppedData');
range = [min(vol(:)) max(vol(:))];

sliderSize = 0.02*diff(range);

% create colorbar for the first time
handles.ColorbarImage = imagesc(linspace(clim(1),clim(2),steps),...
	        0.5,...
			reshape(cmap,1,steps,3));
set(handles.Colorbar,'ylim',[0 1],...
	'ytick',[],...
	'xlim',range,...
	'xDir','normal',...
	'ydir','normal',...
	'FontSize',8,...
	'tag','Colorbar');

guidata(handles.volViewerFigure,handles);




%% SETCOLORBAR function
function ch = setColorbar(handles)
[~, handles] = updateVolViewerColorBar(handles);
ch = handles.Colorbar;
setappdata(ch,'motionpointer','SOM top');