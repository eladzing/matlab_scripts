% ARROW - adds an arrow to an axis.
%
% Usage:
% ------
% p=arrow(parent,dir,pos)
% Inputs are:
%   parent - the object on which to place the arrow
%   dir - the direction of the arrow (down, up, left, right)
%   pos - the position on the axis to place the arrow
% Outputs are:
%   p - pointer to the newly created arrow

% By Ran Klein 1-Nov-2005

% Arrow points are as follows:
%   21012    21012      12345     12345
% 5  *-*   5   *     2   *     2   *
% 4  | |   4  / \    1 *-*\    1  /*-*
% 3 ** **  3 ** **   0 |   *   0 *   |
% 2  \ /   2  | |   -1 *-*/   -1  \*-*
% 1   *    1  *-*   -2   *    -2   *

function p=arrow(parent,dir,pos)


switch dir
	case 'down'
		pts=[ 0 0; -2 0.5; -1 0.5; -1 1; 1 1; 1 0.5; 2 0.5 ];
		mp = 'SOM leftright';
	case 'up'
		pts=[ 0 1; 2 0.5; 1 0.5; 1 0; -1 0; -1 0.5; -2 0.5; ];
		mp = 'SOM leftright';
	case 'right'
		pts=[ 1 0; 0.5 -2; 0.5 -1; 0 -1; 0 1; 0.5 1; 0.5 2 ];
		mp = 'SOM topbottom';
	case 'left'
		pts=[ 0 0; 0.5 2; 0.5 1; 1 1; 1 -1; 0.5 -1; 0.5 -2 ];
		mp = 'SOM topbottom';
end

% Modify the arrows to look good no matter what
% the data aspect ratio may be.
if strcmp(dir,'up') || strcmp(dir,'down')
	lim=get(parent,'xlim');
	scale=abs(lim(1)-lim(2))/15/5;
	pts(:,1)=pts(:,1)*scale+pos;
else
	lim=get(parent,'ylim');
	scale=abs(lim(1)-lim(2))/15/5;
	pts(:,2)=pts(:,2)*scale+pos;
end

% Create the patches, and add app data to them to remember what
% They are associated with.
handles=guidata(parent);
h=getappdata(handles.volViewerFigure,'SliceContextMenu');
p=patch('vertices',pts,'faces',1:size(pts,1),...
	'facec','g','facea',.5,'edgec','k','linewidth',2,...
	'hittest','on','buttondownfcn','objectMove(gcbo)',...
	'parent',parent,'tag','SliceArrow',...
	'UIcontextmenu',h.sliceMenu,...
	'handlevisibility','on');
setappdata(p,'CenterPosition',pos);
setappdata(p,'motionpointer',mp);
