% FRAMESELECTIONCHANGE - handles a change in selection of time frames in
% the volViewer interface.

% By Ran Klein, University of Ottawa Heart Institute, 27-Oct-2005


function FrameSelectionChange(hObject,eventdata,handles)

if nargin<1 
	handles = guidata(gcf); % Get the handle to all the other handles in the TimeFramePanel
else % call back triggered this call then update GUI objects
	if nargin<3
		handles = guidata(hObject);
	end
	switch get(hObject,'Tag') % Based on the triggering object do
		case 'SingleFrame'
			set([handles.EndFrameSlider,handles.EndFrame,handles.FrameString],'Enable','Off');
			set([handles.StartFrameSlider,handles.StartFrame], 'Enable','On');
			set(handles.SingleFrame,'Value',1);
			set([handles.FrameRange handles.Other],'Value',0);
		case 'FrameRange'
			set([handles.FrameString],'Enable','Off');
			set([handles.StartFrameSlider,handles.StartFrame,handles.EndFrameSlider,handles.EndFrame], 'Enable','On');
			set(handles.FrameRange,'Value',1);
			set([handles.SingleFrame handles.Other],'Value',0);
		case 'Other'
			set([handles.StartFrameSlider,handles.StartFrame,handles.EndFrameSlider,handles.EndFrame], 'Enable','Off');
			set([handles.FrameString],'Enable','On');
			set(handles.Other,'Value',1);
			set([handles.SingleFrame handles.FrameRange],'Value',0);
		case 'StartFrameSlider'
			first = round(get(handles.StartFrameSlider,'Value'));
			set(handles.StartFrameSlider,'Value',first);
			set(handles.StartFrame,'String',first);
			last = round(get(handles.EndFrameSlider,'Value'));
			last = max(last, first);
			set(handles.EndFrameSlider,'Value',last);
			set(handles.EndFrame,'String',last);
		case 'EndFrameSlider'
			last = round(get(handles.EndFrameSlider,'Value'));
			set(handles.EndFrameSlider,'Value',last);
			set(handles.EndFrame,'String',last);
			first = round(get(handles.StartFrameSlider,'Value'));
			first = min(last, first);
			set(handles.StartFrameSlider,'Value',first);
			set(handles.StartFrame,'String',first);
		case 'FrameString'
			set([handles.StartFrameSlider,handles.StartFrame,handles.EndFrameSlider,handles.EndFrame], 'Enable','Off');
			set([handles.FrameString],'Enable','On');
		case {'volViewerFigure','TimeOperationPopupmenu'}
			% Do nothing - this was just a call to update the display.
		otherwise
			msgId = sprintf('volViewer:%s:unrecognizedObjectHandle',mfilename);
			error(msgId,['An unexpected object handle (Tag=' get(hObject,'Tag') ') triggered this routine.']);
	end % switch
end


% Recalculate the data to display
% ===============================
set(handles.volViewerFigure,'Pointer','watch'); % Set mouse pointer
drawnow;

data = getappdata(handles.DisplayAxis,'RawData');
if isappdata(handles.volViewerFigure,'TimePoints')
	frame_length = mst2frameTimes(getappdata(handles.volViewerFigure,'TimePoints'));
else
	frame_length = ones(1,size(data,4));
end
if get(handles.SingleFrame,'Value') % Use a single frame
	if length(size(data)) == 3 % Single frame mode
		d = data;
	else
		frame = round(get(handles.StartFrameSlider,'Value'));
		d = data(:,:,:,frame);
		frame_length = frame_length(frame);
	end
elseif get(handles.FrameRange,'Value')
	first = round(get(handles.StartFrameSlider,'Value'));
	last = round(get(handles.EndFrameSlider,'Value'));
	d = data(:,:,:,first:last);
	frame_length = frame_length(first:last);
elseif get(handles.Other,'Value')
	mask = filterString2Mask(get(handles.FrameString,'String'),size(data,4));
	if ischar(mask) % invalid mask
		set(handles.FrameString,'String',mask);
	else
		d = data(:,:,:,mask);
		frame_length = frame_length(mask);	end
else
	return % No frame mode has been selected - ignore
end

switch upper(getPullDownValue(handles.TimeOperationPopupmenu))
	case 'SUM'
		d = sum(double(d),4);
	case 'AVERAGE'
		d = mean(double(d),4);
	case 'WEIGHTED AVERAGE'
		d = sum(scaleFrames(double(d),frame_length),4)/sum(frame_length);
	case 'INTEGRATE'
		d = sum(scaleFrames(double(d),frame_length),4);
	case 'MAX'
		d = max(double(d),[],4);
	otherwise
		error(['Unknown operation: ' getPullDownValue(handles.TimeOperationPopupmenu)]);
end
setappdata(handles.DisplayAxis,'SummedData',double(d)); % Double-d - that's huge!!!!

% Crop the new data and then update the display
cropChange(handles);

set(handles.volViewerFigure,'Pointer','arrow')