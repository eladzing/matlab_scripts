% Handle generic motion events for the figure window.

% By Ran Klein, University of Ottawa Heart Institute, 31-Oct-2005


function volViewerPointerMotion(hObject, event, handles)

handles=guidata(hObject);
% figure(handles.volViewerFigure)
hObject = hittest(handles.volViewerFigure);
hittestAxes = getParentAxes(hObject);

% Some objects get a special pointer when the mouse waves
% over them.  Get it from appdata.
if ~isempty(hObject)
	if ~getappdata(handles.DisplayAxis,'BlockSliceControl')
		mp = getappdata(hObject,'motionpointer'); % Get the pointer pattern for the object
		ph = get(handles.volViewerFigure,'pointer'); % Get the handle to the pointer

		if isempty(mp) % If there is no specific pointer for this object
			mp = get(0,'defaultfigurepointer'); % Use the default
		end

		if isa(mp,'char') && isa(ph,'char') && ~strcmp(mp,ph)
			setpointer(handles.volViewerFigure, mp);  % Set the pointer
		end
	end
end

% Update the current control object - if one exists
cObj = getappdata(handles.DisplayAxis,'ControlObject');
if ~isempty(cObj)
	ph = get(cObj.handle,'Parent'); % Get the parent handle of the current object
	% Slice control popped over from one axis to another?
	if ~isempty(hittestAxes) && ph~=hittestAxes && ...
			strcmp(get(cObj.handle,'Tag'),'SliceArrow') && ph~=handles.Colorbar && ...
			any(hittestAxes==[handles.XAxis, handles.YAxis, handles.ZAxis])
		pos = get(cObj.tipHandle,'Position');
		% rotate/flip arrows and rescale, set center position, fix tooltips, adjust mouse pointer
		switch hittestAxes
			case handles.XAxis
				if ph==handles.YAxis
					set(cObj.handle,'xdata',get(cObj.handle,'ydata')/diff(get(handles.YAxis,'ylim'))*diff(get(handles.XAxis,'xlim')),'ydata',1-get(cObj.handle,'xdata'));
				else % ZAxis
					set(cObj.handle,'xdata',get(cObj.handle,'ydata')/diff(get(handles.ZAxis,'ylim'))*diff(get(handles.XAxis,'xlim')),'ydata',get(cObj.handle,'xdata'));
				end
				mp = 'SOM leftright';
				pos(2) = -0.5;
				set(cObj.tipHandle,'VerticalAlignment','top', 'HorizontalAlignment','center');
				centerPos = mean(get(cObj.handle,'xdata'));
			case handles.YAxis
				if ph==handles.XAxis
					set(cObj.handle,'xdata',1-get(cObj.handle,'ydata'),'ydata',get(cObj.handle,'xdata')/diff(get(handles.XAxis,'xlim'))*diff(get(handles.YAxis,'ylim')));
				else % ZAxis
					set(cObj.handle,'xdata',1-get(cObj.handle,'xdata'),'ydata',get(cObj.handle,'ydata')/diff(get(handles.ZAxis,'ylim'))*diff(get(handles.YAxis,'ylim')));
				end
				mp = 'SOM topbottom';
				pos(1) = 1.5;
				set(cObj.tipHandle,'VerticalAlignment','middle', 'HorizontalAlignment','left');
				centerPos = mean(get(cObj.handle,'ydata'));
			case handles.ZAxis
				if ph==handles.YAxis
					set(cObj.handle,'xdata',1-get(cObj.handle,'xdata'),'ydata',get(cObj.handle,'ydata')/diff(get(handles.YAxis,'ylim'))*diff(get(handles.ZAxis,'ylim')));
				else % XAxis
					set(cObj.handle,'xdata',get(cObj.handle,'ydata'),'ydata',get(cObj.handle,'xdata')/diff(get(handles.XAxis,'xlim'))*diff(get(handles.ZAxis,'ylim')));
				end
				mp = 'SOM topbottom';
				pos(1) = -0.5;
				set(cObj.tipHandle,'VerticalAlignment','middle', 'HorizontalAlignment','right');
				centerPos = mean(get(cObj.handle,'ydata'));
		end
		ph = hittestAxes;
		set(cObj.handle,'parent',ph);
		setpointer(handles.volViewerFigure,mp)
		setappdata(cObj.handle,'CenterPosition',centerPos)
		setappdata(cObj.handle,'motionpointer',mp);
		set(cObj.tipHandle,'Position',pos,'Parent',ph);
		drawnow;
	end
	
	pos=get(ph,'currentpoint'); % Get the position of the point on the active axis

	value=getappdata(cObj.handle,'CenterPosition');	% the line the arrow points at.

	% Bind the axes position to the limits of that axes.
	switch get(cObj.handle,'Tag')
		case 'MinRangeBar'
			sliceh = findobj(get(ph,'Children'),'Tag','SliceArrow');
			limit = zeros(1, length(sliceh));
			for i=1:length(limit)
				limit(i) = getappdata(sliceh(i),'CenterPosition');
			end
			limit = [1 min([getappdata(ph,'MaxRange'); limit])];
		case 'MaxRangeBar'
			sliceh = findobj(get(ph,'Children'),'Tag','SliceArrow');
			limit = zeros(1, length(sliceh));
			for i=1:length(limit)
				limit(i) = getappdata(sliceh(i),'CenterPosition');
			end
			limit = [max([getappdata(ph,'MinRange'); limit]), getappdata(ph,'MaxMaxRange')];
		case 'SliceArrow'
			if ph==handles.Colorbar
				limit = get(ph,'xlim');
			else
				limit = [getappdata(ph,'MinRange'), getappdata(ph,'MaxRange')];
			end
		otherwise
			error('Unidentified control object detected by volViewerPointerMotion');
	end

	if ph==handles.YAxis || ph==handles.ZAxis% Y or Z axes
		pos(1,2) = min(max(round(pos(1,2)),limit(1)),limit(2));
		ydiff=pos(1,2)-value;
		v=get(cObj.handle,'vertices');
		v(:,2)=v(:,2)+ydiff;
		set(cObj.handle,'vertices',v);
		if ph==handles.YAxis % Y Axis
			if pos(1,2)~=round(value)
				updateSlice(handles, [], pos(1,2), [], cObj.slice);
			end
			tippos = [1.5 value];
		else  % ZAxis
			if pos(1,2)~=round(value)
				updateSlice(handles, [], [], pos(1,2), cObj.slice);
			end
			tippos = [-0.5 value];
		end
		setappdata(cObj.handle,'CenterPosition',pos(1,2));
		value = pos(1,2);
	else %if ph==handles.XAxis% Colorbar or the X Axis
		if ph==handles.XAxis
			pos(1,1) = min(max(round(pos(1,1)),limit(1)),limit(2));
		else
			pos(1,1) = min(max(pos(1,1),limit(1)),limit(2));
		end
		xdiff=pos(1,1)-value;
		v=get(cObj.handle,'vertices');
		v(:,1)=v(:,1)+xdiff;
		set(cObj.handle,'vertices',v);
		if  ph==handles.XAxis 
			if pos(1,1)~=round(value)
				updateSlice(handles, pos(1,1), [], [], cObj.slice);
			end
			tippos = [value -0.5];
		elseif ph==handles.Colorbar
			updateContourVals(handles.Colorbar);
			updateContour(handles);
			tippos = [value -0.5];
		end
		setappdata(cObj.handle,'CenterPosition',pos(1,1));
		value = pos(1,1);
	end

	if ph==handles.XAxis || ph==handles.YAxis || ph==handles.ZAxis
		set(cObj.tipHandle,'string',sprintf('Value: %1.0f',round(value)),...
			'position', tippos, ...
			'visible','on');
	else
		set(cObj.tipHandle,'string',sprintf('Value: %1.2g',value),...
			'position', tippos, ...
			'visible','on');
	end
	% Do not show the metaslice
	set(getappdata(handles.DisplayAxis,'metaSlice'),'visible','off');
else % Not dragging an object

	% Pointer is not over one of the axis scales then do not show the meta
	% slice
	if isempty(hObject) || (hObject ~= handles.XAxis && hObject ~= handles.YAxis && hObject ~= handles.ZAxis)
		set(getappdata(handles.DisplayAxis,'metaSlice'),'visible','off');
	else % hObject is hovering over an axis
		if ~getappdata(handles.DisplayAxis,'BlockSliceControl')
			pos = get(hObject,'currentpoint');

			xl = get(handles.DisplayAxis,'xlim'); % Get display coordinate limits
			yl = get(handles.DisplayAxis,'ylim');
			zl = get(handles.DisplayAxis,'zlim');

			switch hObject
				case handles.XAxis
					xdata = [ pos(1,1) pos(1,1) pos(1,1) pos(1,1) pos(1,1) ];
					ydata = [ yl(1) yl(2) yl(2) yl(1) yl(1) ];
					zdata = [ zl(2) zl(2) zl(1) zl(1) zl(2) ];
				case handles.YAxis
					ydata = [ pos(1,2) pos(1,2) pos(1,2) pos(1,2) pos(1,2) ];
					xdata = [ xl(1) xl(2) xl(2) xl(1) xl(1) ];
					zdata = [ zl(2) zl(2) zl(1) zl(1) zl(2) ];
				case handles.ZAxis
					zdata = [ pos(1,2) pos(1,2) pos(1,2) pos(1,2) pos(1,2) ];
					ydata = [ yl(1) yl(2) yl(2) yl(1) yl(1) ];
					xdata = [ xl(2) xl(2) xl(1) xl(1) xl(2) ];
			end

			set(getappdata(handles.DisplayAxis,'metaSlice'),'Parent',handles.DisplayAxis,...
				'xdata',xdata,'ydata',ydata,'zdata',zdata,'visible','on');
			drawnow;
		end
	end
end



%% HELPER FUNCTIONS

%% Get the axes that the object belongs to.
function hObject = getParentAxes(hObject)
if strcmpi(get(hObject,'type'),'figure')
	hObject = [];
else
	if ~strcmpi(get(hObject,'type'),'axes')
		hObject = get(hObject,'parent');
	end
end