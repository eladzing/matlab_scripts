% ControlBar - adds an arrow to an axis.
%
% Usage:
% ------
% p=arrow(parent,dir,pos,tag,color)
% Inputs are:
%   parent - the object on which to place the control
%   dir - the direction of the patch arrow (down, up, left, right)
%   pos - the position on the axis to place the conrol patch
%   tag - the tag name for the control patch
%   color - the color of the control patch
% Outputs are:
%   p - pointer to the newly created arrow

% By Ran Klein, University of Ottawa Heart Institute, 2014-01-15

function p = ControlBar(parent,dir,pos,tag,color)

if nargin<4
	tag = 'ControlBar';
end
if nargin<5
	color = 'g';
end

switch dir
	case 'vertical'
		pts=[ 0  1; 2 0.5; 0 0; -2 0.5];
		mp = 'SOM leftright';
		% Modify the arrows to look good no matter what the data aspect
		% ratio may be.
		lim=get(parent,'xlim');
		scale=abs(lim(1)-lim(2))/40/3;
		pts(:,1)=pts(:,1)*scale+pos;
	case 'horizontal'
		pts=[ 1 0; 0.5 2; 0 0; 0.5 -2];
		mp = 'SOM topbottom';
		% Modify the arrows to look good no matter what the data aspect
		% ratio may be.
		lim=get(parent,'ylim');
		scale=abs(lim(1)-lim(2))/40/3;
		pts(:,2)=pts(:,2)*scale+pos;
end

% Create the patches, and add app data to them to remember what
% They are associated with.
p = findobj(get(parent,'children'),'flat','tag',tag);
if isempty(p)
	p=patch('vertices',pts,'faces',1:size(pts,1),...
		'facec',color,'edgec','k',...
		'hittest','on','buttondownfcn','objectMove(gcbo)',...
		'parent',parent,'tag',tag,...
		'handlevisibility','on');
else
	set(p,'vertices',pts,'faces',1:size(pts,1),...
		'facec',color,'edgec','k',...
		'hittest','on','buttondownfcn','objectMove(gcbo)',...
		'handlevisibility','on');
end
setappdata(p,'CenterPosition',pos);
setappdata(p,'motionpointer',mp);
