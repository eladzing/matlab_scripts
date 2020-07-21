% VOLVIEWER M-file for volViewer.fig
%
% volviewer is a 4D data set viewing application. Three axes control
% viewing slices in the volume and end sliders can be used to crop the
% data. A set of sliders in the top left may be used to control the 4th
% dimension (time/factor slice). Refer to the help option in the volViewer
% menu for additional information.
% 
% volviewer(data) - view the data (a 3D/4D matrix).
%
% volviewer(data, subsetstr) - mask string for the slice/volume to display. 
%      Defines either a subspace in the volume of data to crop the data
%      and/or slices with in the volume that should be displayed at
%      startup. By default, no cropping or slices are applied. Refer to
%      filterSubset2Mask for more details on the string format.
%
% volviewer(data, subsetestr, timestr) - also sets the 4thD time/factors
%      to display. Refer to filterString2Mask for more details on the
%      string format.
%
% volviewer(data, subsetstr, timestr, propertie, value, ...) - Pairs of 
%      property and values specific to the volViewer: Properties include:
%  'FigureName' - String of the name of the figure.
%  'FigureColor' - Colorspec to set the background color of th figure.
%  'Position' - The position on the screen for the window can be specified
%      in two forms:
%      - [x y width height] - in normalized units [0-1].
%      - a string in the format supported by the movegui function.
%  'FramePanelTitle' - String of title of the frame panel. If not
%      specified. Default is 'Time Frames'.
%  'AxisNames' - Cell of strings (3-elemnet) of the names of each of first 
%      3-dimension. Default is {'x-axis','y-axis','z-axis'}.
%  'BlockSliceControl' - Logical value of whether to block or allow 
%      control of the slices to those provided by subspacestr. If not
%      specified default is false.
%  'PixelDimensions' - Dimensions of the pixels [dx dy dz] to control the
%      aspect ratio.
%  'Units' - string name of the image units. Default is empty.
%  'TimeUnits' - string name of time units. Default is empty.
%  'Time' - The middle frame times
%  'TimeOp' - Operation to carry out on time frames
%            [{'Sum','Average','Weighted Average','Integrate','Max']
%  'Colorbar' - Logical value of whether to show the colorbar or not.
%      If not specified, default is true.
%  'Contours' - Array of scalars of the intensity levels for which
%      contour lines are display. If not specified default is contour at
%      zero.
%  'ExtraSurface' - A surface, whose contours will be superimposed on the
%       slices. Extra surface is a structure or a cell of structures with
%       the following form:
%          struct('X',xdata, 'Y',ydata, 'Z',zdata, 'LineStyle','.r', 'Type','plot')
%       where the X, Y and Z data are the coordinates of the surface in 
%       pixel units and can be cells for each time frame or non-cell if 
%       common for all time frames, LineStyle is the color and plot styles 
%       of the plot and Type:
%       - plot - plot individual points at XYZ.
%       - mesh - plot lines of the intersection of the surface with the
%       slices.
%  'Colormap' - sets the colormap to be used. can be either a string to
%      one of the supported colormaps in, or an array definining the
%      colormap. If a string is used, it may be preceded by the letter "i"
%      to indicate inversion on the colormap order.
%  'PointerCallback' - specifies a callback function to execute when
%      user clicks on slices. The function ÿvolViewerCoordÿ may be used to 
%      determine the coordinates of the point of overlap of the mouse 
%      pointer and the selected object. See DisplayPixelCurve.m for an
%      example.
%  'WaitForClose' - If no output parameters are specified, forces the 
%      execution to wait until the volviewer is closed. If output
%      parameters do exist, this field is ignored.
%
% str = volviewer(...) - returns a subset selection string in the format of
%      subspacestr, whos contents depends on the slices, crops, and the
%      acquire button pressed. (Nan if user chooses to not return the value)
% [str, frames] = volviewer(...) - Also returns the frames selected in the
%      GUI. (Nan if user chooses to not return the value)
%
% See also: filterSubset2Mask, filterString2Mask, DisplayPixelCurve,
% volViewerDemo
%
% User manual included: volViewer Manual.pdf
%
% By Ran Klein, University of Ottawa Heart Institute, 2014-01-15
% 2019-02-15 - RK - Updated to support Matlab version 2017b.


function varargout = volViewer(varargin)

% Last Modified by GUIDE v2.5 15-Nov-2011 15:23:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @volViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @volViewer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if length(varargin)<2
	varargin{2} = '';
end
if length(varargin)<3
	varargin{3} = '';
end

fig = findobj(get(0,'children'),'flat','tag','volViewerFigure');
if gui_Singleton && ~isempty(fig) && strcmpi(get(fig,'visible'),'on') && ~ischar(varargin{1})
	initializeVolViewer(guidata(fig), varargin);
else
	if nargout
		varargin = [varargin, 'Nargout', nargout];
		[varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
	else
		gui_mainfcn(gui_State, varargin{:});
	end
end
% End initialization code - DO NOT EDIT

% --- Executes just before volViewer is made visible.
function volViewer_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for volViewer
handles.output = {nan,nan};
guidata(hObject, handles); % Update handles structure

% get height of panels for future resize
units = get(handles.TimeFramePanel,'units');
set(handles.TimeFramePanel,'units','pixel');
pos = get(handles.TimeFramePanel,'pos');
setappdata(handles.TimeFramePanel,'height',pos(4));
set(handles.TimeFramePanel,'units',units);

if ~isempty(varargin) % Keep the display data
	initializeVolViewer(handles, varargin);
else
	error('Data not passed to volViewer')
end
set(handles.PlayButton,'BackgroundColor',[0 0.3 0]);
OneWayLoopButton_Callback(handles.OneWayLoopButton,[],handles);
drawnow;

% Wait until the visibility of the figure is turned off
if getappdata(handles.volViewerFigure,'WaitForClose')
	waitfor(handles.volViewerFigure,'Visible');
	try
		delete(getappdata(gcf,'ChildObjects'));
	catch
	end
end


% --- Outputs from this function are returned to the command line.
function varargout = volViewer_OutputFcn(hObject, eventdata, handles)
varargout = handles.output; 
if getappdata(handles.volViewerFigure,'WaitForClose')
	delete(handles.volViewerFigure); % Close the window once done.
else
	pos = getappdata(handles.volViewerFigure,'DefPos');
	if ~isempty(pos)
		if ischar(pos)
			movegui(handles.volViewerFigure ,pos);
		else
			set(handles.volViewerFigure,'units','normalized','position',pos);
		end
	end
end


% --- Executes on button press in AcquireSlices.
function AcquireSlices_Callback(hObject, eventdata, handles)
h=findobj([handles.XAxis, handles.YAxis, handles.ZAxis],'Type','patch','Tag','SliceArrow');
if ~isempty(h) % Slices are indicated
	coord = getappdata(handles.volViewerFigure,'MaskCoordMode');
	if isempty(coord) % default coordination mode
		coord = {'x','y','z'};
	end
	str=[];
	for i=1:length(h)
		switch get(h(i),'Parent')
			case handles.XAxis
				str = [str ',' coord{1} num2str(getappdata(h(i),'CenterPosition'))];
			case handles.YAxis
				str = [str ',' coord{2} num2str(getappdata(h(i),'CenterPosition'))];
			case handles.ZAxis
				str = [str ',' coord{3} num2str(getappdata(h(i),'CenterPosition'))];
		end
	end
	if length(str)>1 % Remove the first , in the string
		str=str(2:end);
	end
else
	str = [];
end

if get(handles.AcquireTimeFrames,'Value')==1
	% Get the frames that were used
	if get(handles.SingleFrame,'Value') % Single frame
		frames = num2str(round(get(handles.StartFrameSlider,'Value')));
	elseif get(handles.FrameRange,'Value') % Frame Range
		frames = [num2str(round(get(handles.StartFrameSlider,'Value'))) '-'...
			num2str(round(get(handles.EndFrameSlider,'Value')))];
	elseif get(handles.Other,'Value') % Frame string
		frames = get(handles.FrameString,'String');
	end
else
	frames = nan;
end

handles.output = {str, frames};
guidata(hObject, handles);
set(handles.volViewerFigure,'Visible','off');


% --- Executes on button press in AcquireVolume.
function AcquireVolume_Callback(hObject, ~, handles)
% hObject    handle to AcquireVolume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if getappdata(handles.XAxis,'MinRange')~=1 || ...
		getappdata(handles.XAxis,'MaxRange')~=getappdata(handles.XAxis,'MaxMaxRange') || ...
		getappdata(handles.YAxis,'MinRange')~=1 || ...
		getappdata(handles.YAxis,'MaxRange')~=getappdata(handles.YAxis,'MaxMaxRange') || ...
		getappdata(handles.ZAxis,'MinRange')~=1 || ...
		getappdata(handles.ZAxis,'MaxRange')~=getappdata(handles.ZAxis,'MaxMaxRange')

	str = ['(' num2str(getappdata(handles.XAxis,'MinRange')) '-'...
		num2str(getappdata(handles.XAxis,'MaxRange')) ',',...
		num2str(getappdata(handles.YAxis,'MinRange')) '-'...
		num2str(getappdata(handles.YAxis,'MaxRange')) ',',...
		num2str(getappdata(handles.ZAxis,'MinRange')) '-'...
		num2str(getappdata(handles.ZAxis,'MaxRange')) ')'];
else
	str = '';
end


if get(handles.AcquireTimeFrames,'Value')==1
	% Get the frames that were used
	if get(handles.SingleFrame,'Value') % Single frame
		frames = num2str(round(get(handles.StartFrameSlider,'Value')));
	elseif get(handles.FrameRange,'Value') % Frame Range
		frames = [num2str(round(get(handles.StartFrameSlider,'Value'))) '-'...
			num2str(round(get(handles.EndFrameSlider,'Value')))];
	elseif get(handles.Other,'Value') % Frame string
		frames = get(handles.FrameString,'String');
	end
else
	frames = nan;
end

handles.output = {str, frames};
guidata(hObject, handles);
set(handles.volViewerFigure,'Visible','off');


% --- Executes on button press in AcquireBoth.
function AcquireBoth_Callback(hObject, eventdata, handles)
h=findobj([handles.XAxis, handles.YAxis, handles.ZAxis],'Type','patch','Tag','SliceArrow');
if ~isempty(h) % Slices are indicated
	coord = getappdata(handles.volViewerFigure,'MaskCoordMode');
	if isempty(coord) % default coordination mode
		coord = {'x','y','z'};
	end
	str=[];
	for i=1:length(h)
		switch get(h(i),'Parent')
			case handles.XAxis
				str = [str ',' coord{1} num2str(getappdata(h(i),'CenterPosition'))];
				ds = [get(handles.YAxis,'ylim') get(handles.ZAxis,'ylim')];
				cs = [getappdata(handles.YAxis,'MinRange') getappdata(handles.YAxis,'MaxRange'),...
					getappdata(handles.ZAxis,'MinRange') getappdata(handles.ZAxis,'MaxRange')];
				if any(ds~=cs)
					str=sprintf('%s(%d:%d,%d:%d)',str,cs(1),cs(2),cs(3),cs(4));
				end
			case handles.YAxis
				str = [str ',' coord{2} num2str(getappdata(h(i),'CenterPosition'))];
				ds = [get(handles.XAxis,'xlim') get(handles.ZAxis,'ylim')];
				cs = [getappdata(handles.XAxis,'MinRange') getappdata(handles.XAxis,'MaxRange'),...
					getappdata(handles.ZAxis,'MinRange') getappdata(handles.ZAxis,'MaxRange')];
				if any(ds~=cs)
					str=sprintf('%s(%d:%d,%d:%d)',str,cs(1),cs(2),cs(3),cs(4));
				end
			case handles.ZAxis
				str = [str ',' coord{3} num2str(getappdata(h(i),'CenterPosition'))];
				ds = [get(handles.XAxis,'xlim') get(handles.YAxis,'ylim')];
				cs = [getappdata(handles.XAxis,'MinRange') getappdata(handles.XAxis,'MaxRange'),...
					getappdata(handles.YAxis,'MinRange') getappdata(handles.YAxis,'MaxRange')];
				if any(ds~=cs)
					str=sprintf('%s(%d:%d,%d:%d)',str,cs(1),cs(2),cs(3),cs(4));
				end
		end
	end
	if length(str)>1 % Remove the first , in the string
		str=str(2:end);
	end
else % Cropped volume
	if getappdata(handles.XAxis,'MinRange')~=1 || ...
			getappdata(handles.XAxis,'MaxRange')~=getappdata(handles.XAxis,'MaxMaxRange') || ...
	   getappdata(handles.YAxis,'MinRange')~=1 || ...
			getappdata(handles.YAxis,'MaxRange')~=getappdata(handles.YAxis,'MaxMaxRange') || ...
	   getappdata(handles.ZAxis,'MinRange')~=1 || ...
			getappdata(handles.ZAxis,'MaxRange')~=getappdata(handles.ZAxis,'MaxMaxRange')
		
			str = ['(' num2str(getappdata(handles.XAxis,'MinRange')) '-'...
			              num2str(getappdata(handles.XAxis,'MaxRange')) ',',...
			           num2str(getappdata(handles.YAxis,'MinRange')) '-'...
			              num2str(getappdata(handles.YAxis,'MaxRange')) ',',...
			           num2str(getappdata(handles.ZAxis,'MinRange')) '-'...
			              num2str(getappdata(handles.ZAxis,'MaxRange')) ')'];
	else
		str = '';
	end
end

if get(handles.AcquireTimeFrames,'Value')==1
	% Get the frames that were used
	if get(handles.SingleFrame,'Value') % Single frame
		frames = num2str(round(get(handles.StartFrameSlider,'Value')));
	elseif get(handles.FrameRange,'Value') % Frame Range
		frames = [num2str(round(get(handles.StartFrameSlider,'Value'))) '-'...
			num2str(round(get(handles.EndFrameSlider,'Value')))];
	elseif get(handles.Other,'Value') % Frame string
		frames = get(handles.FrameString,'String');
	end
else
	frames = nan;
end

handles.output = {str, frames};
guidata(hObject, handles);
set(handles.volViewerFigure,'Visible','off');

function AcquireTimeFrames_Callback(hObject, eventdata, handles)
% Nothing to do

% --- Executes on selection change in Colormap.
function Colormap_Callback(hObject, eventdata, handles)
str=get(handles.Colormap,'string');
val=lower(str{get(handles.Colormap,'value')});
if strcmp(val,'custom')
	cmapeditor
else
	cmap=feval(val);
	if get(handles.ColormapInverse,'value')
		cmap = flipud(cmap);
	end
	set(handles.volViewerFigure,'colormap',cmap);
end
updateVolViewerColorBar(handles)

% --- Executes during object creation, after setting all properties.
function Colormap_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ColormapInverse.
function ColormapInverse_Callback(hObject, eventdata, handles)
Colormap_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function EndFrameSlider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function FrameString_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function StartFrame_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function StartFrameSlider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes when volViewerFigure is resized.
function volViewerFigure_ResizeFcn(hObject, eventdata, handles)
if ~isempty(handles)
	height = getappdata(handles.TimeFramePanel,'height');
	h  = [handles.TimeFramePanel, handles.ColormapPanel];
	for i = 1:length(h)
		units = get(h(i),'units');
		set(h(i),'units','pixel');
		pos = get(h(i),'pos');
		pos(2) = pos(2) + pos(4)-height;
		pos(4) = height;
		set(h(i),'pos',pos);
		set(h(i),'units',units);
	end
	h  = [handles.AcquireSlices, handles.AcquireBoth];
	height = height*[0.25 0.75];
	for i = 1:length(h)
		units = get(h(i),'units');
		set(h(i),'units','pixel');
		pos = get(h(i),'pos');
		pos(2) = pos(2) + pos(4)-height(i);
		pos(4) = height(i);
		set(h(i),'pos',pos);
		set(h(i),'units',units);
	end
	h = [handles.AcquireVolume, handles.AcquireTimeFrames];
	posref = [pos(2) pos(2)-height(1)];
	height = height(1)*[1 1];
	for i = 1:length(h)
		units = get(h(i),'units');
		set(h(i),'units','pixel');
		pos = get(h(i),'pos');
		pos(2) = posref(i);
		pos(4) = height(i);
		set(h(i),'pos',pos);
		set(h(i),'units',units);
	end
end



% --- Executes on button press in PlayButton.
function PlayButton_Callback(hObject, eventdata, handles)
if ~strcmpi(get(hObject,'type'),'timer')
	th = getappdata(handles.volViewerFigure,'Timer');
	if isequal(get(handles.PlayButton,'BackgroundColor'),[0 0.3 0]) % Start Play
		if isempty(th)
			th = timer('ExecutionMode','FixedRate','TimerFcn',@PlayButton_Callback,'UserData',handles.volViewerFigure);
% 			th = timer('TimerFcn',View4D('PlayButton_Callback',handles.SingleFrame,[],guidata(handles)));
			setappdata(handles.volViewerFigure,'Timer',th);
		end
		set(handles.PlayButton,'String','Stop','BackgroundColor','r');
		set(th,'Period',calcPeriod(get(handles.FrameRateSlider,'Value')));
		start(th);
	else % pause
		if ~isempty(th)
			stop(th);
		end
		set(handles.PlayButton,'String','Play','backgroundColor',[0 0.3 0]);
	end
else % event triggered by the timer
	tic
	th = hObject;
	handles = guidata(get(th,'UserData'));
	frame = get(handles.StartFrameSlider,'Value');
	if isequal(get(handles.OneWayLoopButton,'BackgroundColor'),[1 0 0])
		if frame>=get(handles.StartFrameSlider,'Max')
			frame = 1;
		else
			frame = frame + 1;
		end
	elseif getappdata(handles.TwoWayLoopButton,'Forward') % back and forth play
		if frame>=get(handles.StartFrameSlider,'Max') % end of sequence
			frame = frame - 1;
			setappdata(handles.TwoWayLoopButton,'Forward',0)
		else
			frame = frame + 1;
		end
	else
		if frame<=1 % start of sequence
			frame = frame + 1;
			setappdata(handles.TwoWayLoopButton,'Forward',1)
		else
			frame = frame - 1;
		end
	end
	set(handles.StartFrameSlider,'Value',frame);
	FrameSelectionChange(handles.SingleFrame);
end


% --- Executes on slider movement.
function FrameRateSlider_Callback(hObject, eventdata, handles)
th = getappdata(handles.volViewerFigure,'Timer');
if ~isempty(th) && strcmpi(get(th,'Running'),'on')
	stop(th);
	set(th,'Period',calcPeriod(get(handles.FrameRateSlider,'Value')));
	start(th);
end

% --- Executes during object creation, after setting all properties.
function FrameRateSlider_CreateFcn(hObject, eventdata, handles)
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in OneWayLoopButton.
function OneWayLoopButton_Callback(hObject, eventdata, handles)
set(hObject,'BackgroundColor','r');
set(handles.TwoWayLoopButton,'BackgroundColor',[0.3 0 0]);
drawnow;


% --- Executes on button press in TwoWayLoopButton.
function TwoWayLoopButton_Callback(hObject, eventdata, handles)
set(hObject,'BackgroundColor','r');
set(handles.OneWayLoopButton,'BackgroundColor',[0.3 0 0]);
drawnow;

function val = calcPeriod(val)
val = max(0.01,round(1000*exp(-val))/1000);


function RecordButton_Callback(hObject, eventdata, handles)
playing = ~isequal(get(handles.PlayButton,'BackgroundColor'),[0 0.3 0]);
if playing % stop play back
	PlayButton_Callback(hObject, eventdata, handles)
end
[file, path] = uiputfile({'*.gif','Graphics Interchange Format (*.gif)';...
	'*.avi','Audio Video Interleaved file (*.avi)'},'Capture video to...','*.gif');

if ischar(file)
			
	blockControl = getappdata(handles.volViewerFigure,'BlockSliceControl');
	setappdata(handles.volViewerFigure,'BlockSliceControl',true);
	h = findobj(handles.volViewerFigure,'type','uicontrol');
	set(h,'enable','off')

	if isequal(get(handles.OneWayLoopButton,'BackgroundColor'),[1 0 0]) % one way loop
		frames = 1:get(handles.StartFrameSlider,'Max');
	else % sweep loop
		frames = [1:get(handles.StartFrameSlider,'Max') get(handles.StartFrameSlider,'Max')-1:-1:2];
	end

	if strcmpi(file(end-3:end),'.avi')
		aviobj = avifile([path filesep file],'FPS',1/calcPeriod(get(handles.FrameRateSlider,'value')),'Quality',95,'VideoName','Video capture from 4D Viewer.');
	elseif strcmpi(file(end-3:end),'.gif')
		cim = [];
	end
	
	for frame = frames
		set(handles.StartFrameSlider,'Value',frame);
		if ~ishandle(handles.volViewerFigure) || strcmpi(get(handles.volViewerFigure,'BeingDeleted'),'on')
			break;
		end
		try % might have error if figure is closed
			FrameSelectionChange(handles.SingleFrame);
			drawnow;
			
			if strcmpi(file(end-3:end),'.avi')
				aviobj = addframe(aviobj,getframe(handles.volViewerFigure));
			elseif strcmpi(file(end-3:end),'.gif')
				f = getframe(handles.volViewerFigure);
				if isempty(cim)
					cim = frame2im(f);
					cim(1,1,1,length(frames))=0;
					fi = 1;
				else
					fi = fi+1;
					cim(:,:,:,fi)=frame2im(f);
				end
			end
		catch
		end
	end
	if strcmpi(file(end-3:end),'.avi')
		aviobj = close(aviobj);
	elseif strcmpi(file(end-3:end),'.gif')
		cim2 = permute(cim,[1 2 4 3]);
		[cim2, cmap] = rgb2ind(reshape(cim2,size(cim,1),[],3),256,'nodither');
		cim2 = reshape(cim2,size(cim,1),size(cim,2),1,size(cim,4));
		imwrite(cim2,cmap,[path filesep file],'gif','DelayTime',calcPeriod(get(handles.FrameRateSlider,'value')),'LoopCount',inf);
	end
	set(h,'enable','on')
	setappdata(handles.volViewerFigure,'BlockSliceControl',blockControl);
end

if playing % restart play back
	PlayButton_Callback(hObject, eventdata, handles)
end



function volViewerFigure_DeleteFcn(hObject, eventdata, handles)
th = getappdata(handles.volViewerFigure,'Timer');
if ~isempty(th) && isvalid(th)
	try
		stop(th);
		delete(th);
	catch
	end
end

% --- Executes on key press with focus on View4DFigure and no controls selected.
function volViewerFigure_KeyPressFcn(hObject, eventdata, handles)
key = get(handles.volViewerFigure,'CurrentCharacter');
switch key
	case char(28) % left arrow
		set(handles.StartFrameSlider,'value', max(1,get(handles.StartFrameSlider,'value')-1) );
		FrameSelectionChange(handles.StartFrameSlider);
	case char(29) % right arrow
		set(handles.StartFrameSlider,'value', min(get(handles.StartFrameSlider,'Max'), get(handles.StartFrameSlider,'value')+1) );
		FrameSelectionChange(handles.StartFrameSlider);
	case char(31) % down arrow
		if ~get(handles.FrameRange,'value')
			FrameSelectionChange(handles.FrameRange);
		else
			set(handles.EndFrameSlider,'value', max(1,get(handles.EndFrameSlider,'value')-1) );
			FrameSelectionChange(handles.EndFrameSlider);
		end
	case char(30) % up arrow
		if ~get(handles.FrameRange,'value')
			FrameSelectionChange(handles.FrameRange);
		else
			set(handles.EndFrameSlider,'value', min(get(handles.EndFrameSlider,'Max'), get(handles.EndFrameSlider,'value')+1) );
			FrameSelectionChange(handles.EndFrameSlider);
		end
end




% --- Executes on selection change in TimeOperationPopupmenu.
function TimeOperationPopupmenu_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function TimeOperationPopupmenu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
