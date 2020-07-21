% VOLVIEWERCOORD - returns the coordinates that the pointer overlaps
% object, hObject.
%
% [x,y,z] = volViewerCoord(hObject) - x, y, and z are the coordinates of
% the point where mouse pointer overlaps with the nearest object in the
% display axis of volViewer.
% [x,y,z,data] = volViewerCoord(hObject) - Also returns the data associated
% with the point at the x,y,z coordinates. The data is an n dimension array
% representing the 4th dimension of thedata.
%
% This function can be used as part of a callback (set with the 
% PopinterCallback property) from VolViewer to determine the coordinates of
% the point that the user has clicked with the mouse.
%
% See also: volViewer

% By Ran Klein, University of Ottawa Heart Institute, 11/9/2006

function [x,y,z,data] = volViewerCoord(hObject)


if strcmp(get(hObject,'Tag'),'Contour') % is the object the contour or the slice?
	hObject = getappdata(hObject,'Slice');
end

ph = get(hObject,'Parent');
pos = round(get(ph,'currentPoint'));
ah = getappdata(hObject,'Arrow');
switch getappdata(hObject,'SliceType')
	case 'X'
		x = getappdata(ah,'CenterPosition');
		y = interp1(pos(:,1),pos(:,2),x);
		z = interp1(pos(:,1),pos(:,3),x);
	case 'Y'
		y = getappdata(ah,'CenterPosition');
		x = interp1(pos(:,2),pos(:,1),y);
		z = interp1(pos(:,2),pos(:,3),y);
	case 'Z'
		z = getappdata(ah,'CenterPosition');
		x = interp1(pos(:,3),pos(:,1),z);
		y = interp1(pos(:,3),pos(:,2),z);
end

if nargout>3
	data = getappdata(ph,'RawData');
	x = min(max(round(x),1),size(data,1));
	y = min(max(round(y),1),size(data,2));
	z = min(max(round(z),1),size(data,3));
	data = reshape(data(x,y,z,:),1,size(data,4));
end