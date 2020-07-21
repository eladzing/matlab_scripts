% cropChange - Update the display axis based on the settings in the GUI 
% objects.
%
% Usage:
% ------
% cropChange(handles) - handles is the handles structure associated with
% the VolumeViewer GUI.

% By Ran Klein, University of Ottawa Heart Institute, 2014-01-15


function cropChange(handles)

if nargin<1 % handles not provided then fetch the children of the active window
	handles = guidata(gcf);
end

d = getappdata(handles.DisplayAxis,'SummedData');

d=d(getappdata(handles.XAxis,'MinRange'):getappdata(handles.XAxis,'MaxRange'),...
	getappdata(handles.YAxis,'MinRange'):getappdata(handles.YAxis,'MaxRange'),...
	getappdata(handles.ZAxis,'MinRange'):getappdata(handles.ZAxis,'MaxRange'));

% Get cropping points
setappdata(handles.DisplayAxis,'SummedCroppedData',d);

set(handles.DisplayAxis,...
	'Xlim', [getappdata(handles.XAxis,'MinRange'), max(1.1,getappdata(handles.XAxis,'MaxRange'))],...
	'Ylim', [getappdata(handles.YAxis,'MinRange'), max(1.1,getappdata(handles.YAxis,'MaxRange'))],...
	'Zlim', [getappdata(handles.ZAxis,'MinRange'), max(1.1,getappdata(handles.ZAxis,'MaxRange'))]);

% The color limits should be updated to optimize contrast w/o saturation
lim=[nanmin(d(:)) nanmax(d(:))];
if lim(1)==lim(2), lim(1) = 0; end
if lim(1)==lim(2), lim(2) = 1; end
set(handles.DisplayAxis,'clim',lim, 'alim',lim);

updateVolViewerColorBar(handles); % update the colorbar too
