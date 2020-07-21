% OBJECTMOVE - sets the object with handel hObject for drag movement with 
% the mouse pointer. This function is called by the buttonDownFcn call back
% of the object that has been clicked. In addition, this function is called
% by the WindowButtonUpFcn when the mouse button is released.

% By Ran Klein, University of Ottawa Heart Institute, 2014-01-14


function objectMove(hObject, event, handles)

if nargin<3
	handles = guidata(gcbf);
end

% Run this function only if the left mouse button was clicked - a right
% click indicates context menu, and there is no need to go through with all
% this.
if ~strcmpi(get(handles.volViewerFigure,'SelectionType'),'normal')
	return
end

% refresh callback to overcome an issue with rotate, zoom, and pan, where
% the callback is overwritten and not restored. This bug can be addressed
% in newer versions of Matlab, but no in version 7 and below.
set(handles.volViewerFigure,'WindowButtonMotionFcn',@volViewerPointerMotion,'WindowButtonUpFcn','objectMove');
cObj = getappdata(handles.DisplayAxis,'ControlObject');
if nargin % A parameters is passed - indicating that an object has been selected
	if ~isempty(cObj) %Is a control object already present?
		rmappdata(handles.DisplayAxis,'ControlObject');
	end
	ph = get(hObject,'Parent');
		
	if ~getappdata(handles.DisplayAxis,'BlockSliceControl') || ph==handles.Colorbar
		cObj.handle = hObject;
		cObj.tipHandle = getappdata(handles.DisplayAxis,'Tip'); % Get the handle to the tip object

		cObj.slice = getappdata(cObj.handle,'Slice'); % retreive the corresponding slice to the object
		% (if none exists, [] is returned resulting in creation of a new slice)

		% Stop playback during crop cahnges
		if contains(get(cObj.handle,'Tag'),'RangeBar') && ~isequal(get(handles.PlayButton,'BackgroundColor'),[0 0.3 0])
			volViewer('PlayButton_Callback',handles.PlayButton,[],handles);
		end
		
		if ph==handles.XAxis % The X Axis
			set(cObj.tipHandle,'verticalalign','top',...
				'horizontalalign','center',...
				'Parent',ph);
			if isempty(cObj.slice)
				cObj.slice = updateSlice(handles, getappdata(cObj.handle,'CenterPosition'), [], [], cObj.slice); % No slice created yet
			end
		elseif ph==handles.YAxis % Y Axis
			set(cObj.tipHandle,'verticalalign','middle',...
				'horizontalalign','left',...
				'Parent',ph);
			if isempty(cObj.slice)
				cObj.slice = updateSlice(handles, [], getappdata(cObj.handle,'CenterPosition'), [], cObj.slice); % No slice created yet
			end
		elseif ph==handles.ZAxis % Z Axis
			set(cObj.tipHandle,'verticalalign','middle',...
				'horizontalalign','right',...
				'Parent',ph);
			if isempty(cObj.slice)
				cObj.slice = updateSlice(handles, [], [], getappdata(cObj.handle,'CenterPosition'), cObj.slice); % No slice created yet
			end
		else % Colorbar
			set(cObj.tipHandle,'verticalalign','top',...
				'horizontalalign','center',...
				'Parent',ph);
			surfh = getappdata(cObj.handle,'Surface');
			if ~isempty(surfh)
				cObj.Surface = true;
				cObj.alpha = get(surfh,'FaceAlpha');
				delete(surfh);
				setappdata(cObj.handle,'Surface',[]);
			else
				cObj.Surface = false;
			end
		end
		cObj.Alpha = get(cObj.slice,'Facealpha');
		set(cObj.slice,'Facealpha',1); % No tranperacy on active control's slice
		setappdata(handles.DisplayAxis,'ControlObject',cObj); % Add reference to the active control object
	end
else
	if ~isempty(cObj)
		rmappdata(handles.DisplayAxis,'ControlObject'); % There is no longer a control object selected
		set(cObj.tipHandle,'Visible','Off','Parent',handles.DisplayAxis); % Remove tip (pointer value) from the screen
		value = getappdata(cObj.handle,'CenterPosition'); % Get the final value
		ph = get(cObj.handle,'Parent');
		if contains(get(cObj.handle,'Tag'),'RangeBar') % If a range bar was moved
			if contains(get(cObj.handle,'Tag'),'Min') % MinRangeBar
				value = max(min(round(value),getappdata(ph,'MaxRange')),1); % Make sure final value is within data size boundaries
				setappdata(ph,'MinRange',value);
			else % MaxRangeBar
				value = max(min(round(value),getappdata(ph,'MaxMaxRange')),getappdata(ph,'MinRange')); % Make sure final value is within data size boundaries
				setappdata(ph,'MaxRange',value);
			end
			setappdata(cObj.handle,'CenterPosition',value);
			delete(getappdata(cObj.slice,'Contour'));
			delete(cObj.slice); % Remove the slice from display
			cropChange;                         % update the cropped data
		else % The control object was a slice indicator
			if ph==handles.Colorbar
				lim = get(ph,'XLim');
				value = max(min(value,lim(2)),lim(1));
				if cObj.Surface
					updateSurface(cObj.handle, handles, cObj.alpha);
				end
			else
				value = max(min(round(value),getappdata(ph,'MaxRange')),getappdata(ph,'MinRange')); % Make sure final value is within data size boundaries
			end
			setappdata(cObj.handle,'CenterPosition',value);
			set(cObj.slice,'Facealpha',cObj.Alpha);
		end
	end
end

if exist('ph','var') && ph==handles.Colorbar
	updateContourVals(handles.Colorbar);
end