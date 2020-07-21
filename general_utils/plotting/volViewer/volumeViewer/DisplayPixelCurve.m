% DisplayPixelCurve - Show the time-activity curve for the pixel that was
% clicked on in volViewer. Example code to generate callback functions.
%
% See also: volViewer, volViewerCoord

% By Ran Klein, University of Ottawa Heart Institute, 2014-01-15


function DisplayPixelCurve(hObject,event,handles)
if nargin<3
	handles = guidata(hObject);
end
volViewerFig = handles.volViewerFigure;

% Creat figure if need be
fig = findobj(get(0,'children'),'tag','PixelCurveDisplay');
if isempty(fig)
	units = get(volViewerFig,'Units');
	set(volViewerFig,'Units','Normalized');
	pos = get(volViewerFig,'Position');
	set(volViewerFig,'Units',units);
	space = [pos(1) pos(2) 1-pos(1)-pos(3) 1-pos(2)-pos(4)];
	switch find(space==max(space),1,'last')
		case 1 % position to left
			pos(3) = max(0.2,pos(1));
			pos(1) = 0;
		case 2 % position below
			pos(4) = max(0.2, pos(2));
			pos(2) = 0;
		case 3 % position to right
			pos(1) = pos(1)+pos(3);
			pos(3) = max(0.2,1-pos(1));
		case 4 % position above
			pos(2) = pos(2)+pos(4);
			pos(4) = max(0.2, 1-pos(2));
	end
	fig = figure;
	drawnow;
	set(fig,'tag','PixelCurveDisplay',...
		'Units','Normalized','Position',pos);
	setappdata(volViewerFig,'ChildObjects',[fig	getappdata(volViewerFig,'ChildObjects')]);
end
[x,y,z,pCurve] = volViewerCoord(hObject);

% Make the plot
clf(fig); ah = axes; set(ah,'parent',fig);
plot(ah,pCurve); 
title(ah,['(' num2str(x) ', ',num2str(y) ', ' num2str(z) ')']);
xlabel(ah,getappdata(volViewerFig,'TimeUnits'));
ylabel(ah,getappdata(volViewerFig,'Units'));

% Place figures on top
figure(fig);
figure(volViewerFig);

