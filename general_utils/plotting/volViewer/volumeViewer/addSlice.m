% ADDSLICE - This function is called by the callback function of the axes
% and adds a new slice to the display.

% By Ran Klein 28-Oct-2005


function addSlice

h = gcbo;
handles = guidata(h);
p=[];
if ~getappdata(handles.DisplayAxis,'BlockSliceControl')
	pos = round(get(h,'currentpoint'));

	switch get(h,'Tag')
		case 'XAxis'
			if pos(1,1)>=getappdata(h,'MinRange') && pos(1,1)<=getappdata(h,'MaxRange')
				p=arrow(h,'down',pos(1,1));
				s=updateSlice(handles, pos(1,1), [], []);
				setappdata(p,'Slice',s); % Link arrow and slice
				setappdata(s,'Arrow',p);
			end
		case 'YAxis'
			if pos(1,2)>=getappdata(h,'MinRange') && pos(1,2)<=getappdata(h,'MaxRange')
				p=arrow(h,'right',pos(1,2));
				s=updateSlice(handles, [], pos(1,2), []);
				setappdata(p,'Slice',s); % Link arrow and slice
				setappdata(s,'Arrow',p);
			end
		case 'ZAxis'
			if pos(1,2)>=getappdata(h,'MinRange') && pos(1,2)<=getappdata(h,'MaxRange')
				p=arrow(h,'left',pos(1,2));
				s=updateSlice(handles, [], [], pos(1,2));
				setappdata(p,'Slice',s); % Link arrow and slice
				setappdata(s,'Arrow',p);
			end
	end

	if ~isempty(p) % was the slice actually created - in range???
		if ~strcmpi(getappdata(handles.volViewerFigure,'SelectionMode'),'Volume')
			set(handles.AcquireSlices,'Enable',get(handles.AcquireVolume,'Enable'));
		end
		objectMove(p);
	end
end

%% Contour line on colorbar
% contours can be added even if slice control is not enabled
if strcmp(get(h,'Tag'),'Colorbar') % contour line
	pos = get(h,'currentpoint');
	p=arrow(h,'up',pos(1,1));
	setappdata(p,'Surface',[]);
	objectMove(p);
	updateContour(handles);
end

