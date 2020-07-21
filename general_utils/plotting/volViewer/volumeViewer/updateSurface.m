% UPDATESURFACE - updates surface display in volViewer after the contour
% control % has been changed.
%
% Usage:
% ------
% ph = updateSurface(ah, handles, falpha)
% Inputs are:
%    - ah - arrow handle - the handle to the arrow contro surface
%    - handles - the VolvumeViewer object handles structure. (optional)
%    - falpha - face alpha property. (default is 1 - no tranparency)
% Outputs are:
%    - ph - patch object handles associated with the control objects, ah.
%
% See also: objectMove, volViewrMenus

% By Ran Klein , University of Ottawa Heart Institute, 2014-01-14


function ph = updateSurface(ah, handles, falpha)

if nargin<2
	handles = guidata(ah);
end
ph = getappdata(ah,'Surface');

if nargin<3
	if ~isempty(ph)
		falpha = get(ph,'FaceAlpha');
	else
		falpha = 1;
	end
end

if ~isempty(ph)
	delete(ph);
end

[x, y, z] = meshgrid(getappdata(handles.YAxis,'MinRange'):getappdata(handles.YAxis,'MaxRange'),...
	getappdata(handles.XAxis,'MinRange'):getappdata(handles.XAxis,'MaxRange'),...
	getappdata(handles.ZAxis,'MinRange'):getappdata(handles.ZAxis,'MaxRange'));
ph = patch(isosurface(y, x, z, getappdata(handles.DisplayAxis,'SummedCroppedData'), getappdata(ah,'CenterPosition')));
set(ph,'Parent',handles.DisplayAxis,'Tag','Surface');
% isonormals(y,x,z,getappdata(handles.DisplayAxis,'SummedCroppedData'), p)
color = get(handles.volViewerFigure,'Colormap');
color = interp1(1:size(color,1),color,interp1(get(handles.Colorbar,'xlim'),[1, size(color,1)],getappdata(ah,'CenterPosition')));
set(ph, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', falpha);
material(ph,'metal');
h = getappdata(handles.volViewerFigure,'SliceContextMenu');
set(ph,'UIContextMenu',h.sliceMenu);
setappdata(ah,'Surface',ph);
setappdata(ph,'Arrow',ah);