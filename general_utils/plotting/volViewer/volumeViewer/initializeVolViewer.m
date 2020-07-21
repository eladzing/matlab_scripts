% Initialize the volViewer display at the creation of the GUI.
% Inputs are:
%   handles - the handle strucutre to the GUI components.
%   p - the parameters passed by the user, as follows:
%   p{1} - The data to display - as a 3D/4D matrix.
%   p{2} - mask string for the slice/volume to display. If p{2} defines a
%          subspace in the volume of p{1}, this defines cropping and
%          slicing for the display. If not specified no cropping or slicing
%          is enfored.
%   p{3} - Time frames mask string to be summed for the display by default.
%          If not specified, last frame is default.
%   p{4...} - Pairs of property and value: Properties include
%      'FramePanelTitle' - String of title of the frame panel. If not
%      specified, default is 'Time Frames'.
%      'BlockSliceControl' - Logical value of weather to block control of
%      the slices to those provided by p{2}, or allow control. If not
%      specified default is false.

% By Ran Klein, University of Ottawa Heart Institute, 2014-01-14


function initializeVolViewer(handles, p)


narginchk(2,2);

%% Figure
% ======
% Set the renderer to OpenGL if possible.
try
	set(handles.volViewerFigure,'Renderer','OpenGL')
catch
	warning('Failed to find OpenGL renderer.');
end
if ~strcmpi(get(handles.volViewerFigure,'renderer'),'OpenGL')
	warning('Renderer is not OpenGL. Some graphic features, such as transperacy, may not be displayed properly.');
end

%% Window event callbacks
% ======================
wbh = waitbar(0,'Creating volViewer window. Please wait','WindowStyle','Modal');
delete(findobj(get(handles.volViewerFigure,'Children'),'type','uitoolbar'));
cameratoolbar('init'); % initilize the cameratoolbar so that window callbacks are recorded properly
volViewerMenus(handles); % Initialize the menus
set(handles.DisplayAxis,'Tag','DisplayAxis'); % Overcome tag name bug in GUIDE
set(handles.AcquireSlices,'Enable','Off');

%% Store data and setup DisplayAxis
% ================================
% Check and condition data
if islogical(p{1})
	p{1} = double(p{1});
end

% Update gui variables
setappdata(handles.DisplayAxis,'RawData',p{1});
setappdata(handles.DisplayAxis,'SummedData',p{1}(:,:,:,end));
setappdata(handles.DisplayAxis,'SummedCroppedData',p{1}(:,:,:,end));

s = size(p{1});
lim=[nanmin(p{1}(:)) nanmax(p{1}(:))];
if lim(1)==lim(2), lim(1) = 0; end
if lim(1)==lim(2), lim(2) = 1; end

set(handles.DisplayAxis,'box','on',...
          'xlim',[1 max(1.1,s(1))],'XDir','Normal',...
		  'ylim',[1 max(1.1,s(2))],'YDir','Normal',...
          'zlim',[1 max(1.1,s(3))],'ZDir','Reverse',...
          'clim',lim,...
          'alim',lim,...
		  'Visible','On');
% axes(handles.DisplayAxis); view(3);
set(handles.DisplayAxis,'View',[-37.5 30]);
setappdata(handles.DisplayAxis,'PixelDimensions',[1 1 1]);
daspect(handles.DisplayAxis,[1 1 1]);
axis(handles.DisplayAxis, 'tight','vis3d');
hold(handles.DisplayAxis,'on');
grid(handles.DisplayAxis,'on');
xlabel(handles.DisplayAxis,'X'); ylabel(handles.DisplayAxis,'Y'); zlabel(handles.DisplayAxis,'Z');

% Create the meta slice line on the main figure
metaSlice = findobj(get(handles.DisplayAxis,'children'),'flat','tag','MetaSlice');
if isempty(metaSlice)
	metaSlice = line('parent',handles.DisplayAxis,...
		'vis','off',...
		'linestyle','--',...
		'marker','none',...
		'linewidth',2,...
		'clipping','off',...
		'Tag','MetaSlice');
end
setappdata(handles.DisplayAxis,'metaSlice',metaSlice);
light('parent',handles.DisplayAxis); lighting(handles.DisplayAxis,'GOURAUD');

%% Time Frames
% ===========
waitbar(0.2,wbh,'Applying frame settings. Please wait');
if length(s)==4 % Is this a time sequence of volumes
	step = 1/(s(4)-1);
	set(handles.StartFrameSlider,'Min',1,'Max',s(4),'Value',1,'SliderStep',[step max(0.1, step)]);
	set(handles.EndFrameSlider,'Min',1,'Max',s(4),'Value',1,'SliderStep',[step max(0.1, step)]);
	if length(p)>2 && ~isempty(p{3}) % Time frames to display have been passed as a parameter
		if ischar(p{3}) % If a string was passed
			timeFrames = filterString2Mask(p{3},s(4)); % convert to mask
		else
			timeFrames = p{3};
		end
		if ischar(timeFrames) % Still a string??? 
			timeFrames = s(4); % Bad string was passed - use last frame
			warning('%s is not a valid filter string. Last frame used in place.' ,p{3});
		elseif length(timeFrames)>1 %Otherwis if an array
			if ~all(timeFrames(2:end)-timeFrames(1:end-1) == 1) % if not a continous range
				timeFrames = p{3}; % then use the string
			end
		elseif isempty(timeFrames)
			timeFrames = s(4);
		end
	else
		timeFrames = s(4); % Default is the last time frames
	end
	
	if ischar(timeFrames) % Other string
		set(handles.SingleFrame,'Value',0);
		set(handles.FrameRange,'Value',0);
		set(handles.Other,'Value',1);
		set(handles.StartFrameSlider,'Value',1,'Enable','Off');
		set(handles.StartFrame,'String',1);
		set(handles.EndFrameSlider,'Value',1,'Enable','Off');
		set(handles.EndFrame,'String',1);
		set(handles.FrameString,'String',timeFrames);
	elseif length(timeFrames)==1 % Single Frame
		set(handles.SingleFrame,'Value',1);
		set(handles.FrameRange,'Value',0);
		set(handles.Other,'Value',0);
		set(handles.StartFrameSlider,'Value',timeFrames,'Enable','On');
		set(handles.StartFrame,'String',timeFrames);
		set(handles.EndFrameSlider,'Value',timeFrames,'Enable','Off');
		set(handles.EndFrame,'String',timeFrames,'Enable','Off');
		set(handles.FrameString,'String','','Enable','Off');
	else  % Range
		set(handles.SingleFrame,'Value',0);
		set(handles.FrameRange,'Value',1);
		set(handles.Other,'Value',0);
		set(handles.StartFrameSlider,'Value',timeFrames(1),'Enable','On');
		set(handles.StartFrame,'String',timeFrames(1));
		set(handles.EndFrameSlider,'Value',timeFrames(end),'Enable','On');
		set(handles.EndFrame,'String',timeFrames(end));
		set(handles.FrameString,'String','','Enable','Off');
	end

	% Callbacks
	set([handles.SingleFrame, handles.FrameRange, handles.Other,...
		handles.StartFrameSlider, handles.EndFrameSlider,...
		handles.FrameString, handles.TimeOperationPopupmenu],...
		'Callback','FrameSelectionChange(gcbo)');
else % Not a time sequence - just a volume image -  no frame selection controls
	set(handles.SingleFrame,'Value',1);
	set(handles.FrameRange,'Value',0);
	set(handles.Other,'Value',0);
	set([handles.StartFrameSlider handles.EndFrame],'Value',0);
	set([handles.StartFrame handles.EndFrame],'String','');
	
	set([handles.SingleFrame, handles.FrameRange, handles.Other,...
		handles.StartFrameSlider, handles.EndFrameSlider,...
		handles.FrameString],...
		'Enable','off');
end

%% Axes Bars and Range Limits
% ==========================
waitbar(0.3,wbh,'Applying crop and slices. Please wait');
set(handles.XAxis,'box','on',...
                  'ytick',[],'xgrid','on','xaxislocation','top',...
                  'xlim',[1 s(1)],'ylim',[0 1],...
                  'layer','top','color','none',...
				  'buttondownfcn','addSlice',...
				  'handlevisibility','off',...
				  'Visible','On');
setappdata(handles.XAxis,'motionpointer','SOM bottom');
setappdata(handles.XAxis,'MaxMaxRange',s(1));

set(handles.YAxis,'box','on',...
                  'xtick',[],'ygrid','on',...
                  'xlim',[0 1],'ylim',[1 s(2)],...
                  'layer','top','color','none',...
				  'buttondownfcn','addSlice',...
				  'handlevisibility','off',...
				  'Visible','On');
setappdata(handles.YAxis,'motionpointer','SOM right');
setappdata(handles.YAxis,'MaxMaxRange',s(2));

set(handles.ZAxis,'box','on',...
                  'xtick',[],'ygrid','on','yaxislocation','right',...
                  'xlim',[0 1],'ylim',[1 max(1.1, s(3))],...
				  'YDir','reverse',...
                  'layer','top','color','none',...
				  'buttondownfcn','addSlice',...
				  'handlevisibility','off',...
				  'Visible','On');
setappdata(handles.ZAxis,'motionpointer','SOM left');
setappdata(handles.ZAxis,'MaxMaxRange',s(3));

% add crop sliders
xlim=[1 s(1)];  ylim=[1 s(2)]; zlim=[1 s(3)];

setappdata(handles.XAxis,'MaxRange',xlim(2));
setappdata(handles.XAxis,'MinRange',xlim(1));
setappdata(handles.YAxis,'MaxRange',ylim(2));
setappdata(handles.YAxis,'MinRange',ylim(1));
setappdata(handles.ZAxis,'MaxRange',zlim(2));
setappdata(handles.ZAxis,'MinRange',zlim(1));

ControlBar(handles.XAxis, 'vertical', xlim(1), 'MinRangeBar', 'g');
ControlBar(handles.XAxis, 'vertical', xlim(2), 'MaxRangeBar', 'r');
ControlBar(handles.YAxis, 'horizontal', ylim(1),  'MinRangeBar','g');
ControlBar(handles.YAxis, 'horizontal', ylim(2), 'MaxRangeBar', 'r');
ControlBar(handles.ZAxis, 'horizontal', zlim(1), 'MinRangeBar', 'g');
ControlBar(handles.ZAxis, 'horizontal', zlim(2), 'MaxRangeBar', 'r');

% Handle predefined slices
if length(p)>1 && ~isempty(p{2})
	mask = filterSubset2Mask(p{2}, s(1:3));
	if ~isempty(mask) && ~ischar(mask)
		% Set the coordinate mode based on the first plane's direction
		switch upper(p{2}(1))
			case {'X','Y','Z'}
				setappdata(handles.volViewerFigure,'MaskCoordMode',{'X' 'Y' 'Z'});
			case {'S','C','T'} 
				if upper(p{2}(1:2)) == 'SH'
					setappdata(handles.volViewerFigure,'MaskCoordMode',{'H' 'V' 'Sh'});
				else
					setappdata(handles.volViewerFigure,'MaskCoordMode',{'S' 'C' 'T'});
				end
			case {'H','V'}, setappdata(handles.volViewerFigure,'MaskCoordMode',{'H' 'V' 'Sh'});
		end
		
		xlim=[s(1) 1];  ylim=[s(2) 1]; zlim=[s(3) 1];
		for i = 1:size(mask,1)
			ph = [];
			if length(mask{i,1})==1 % XAxis slice
				ph=arrow(handles.XAxis,'down',mask{i,1});
				sh=updateSlice(handles, mask{i,1}, [], []);
				ylim(1)=min(ylim(1), mask{i,2}(1));
				ylim(2)=max(ylim(2), mask{i,2}(end));
				zlim(1)=min(zlim(1), mask{i,3}(1));
				zlim(2)=max(zlim(2), mask{i,3}(end));
			elseif length(mask{i,2})==1 % YAxis slice
				ph=arrow(handles.YAxis,'right',mask{i,2});
				sh=updateSlice(handles, [], mask{i,2}, []);
				xlim(1)=min(xlim(1), mask{i,1}(1));
				xlim(2)=max(xlim(2), mask{i,1}(end));
				zlim(1)=min(zlim(1), mask{i,3}(1));
				zlim(2)=max(zlim(2), mask{i,3}(end));
			elseif length(mask{i,3})==1 % ZAxis slice
				ph=arrow(handles.ZAxis,'left',mask{i,3});
				sh=updateSlice(handles, [], [], mask{i,3});
				xlim(1)=min(xlim(1), mask{i,1}(1));
				xlim(2)=max(xlim(2), mask{i,1}(end));
				ylim(1)=min(ylim(1), mask{i,2}(1));
				ylim(2)=max(ylim(2), mask{i,2}(end));
			else % Cropped volume
				xlim(1)=min(xlim(1), mask{i,1}(1));
				xlim(2)=max(xlim(2), mask{i,1}(end));
				ylim(1)=min(ylim(1), mask{i,2}(1));
				ylim(2)=max(ylim(2), mask{i,2}(end));
				zlim(1)=min(zlim(1), mask{i,3}(1));
				zlim(2)=max(zlim(2), mask{i,3}(end));
			end
			if ~isempty(ph)
				setappdata(ph,'Slice',sh); % Link arrow and slice
				setappdata(sh,'Arrow',ph);
				set(handles.AcquireSlices,'Enable','On');
			end
		end % svi loop
% 		xlim=sort(xlim); ylim=sort(ylim); zlim=sort(zlim);
	else
		warning('Subset passed to volViewer is not valid');
	end
end



%% Other options in the function call
% ==================================
% Has a title for the 4th dimension been provided - 'Time Frames' is the
% default
contours = [];
colorbarStatus = 'on';
%% Default parameter values
setappdata(handles.DisplayAxis,'BlockSliceControl',false);
if length(s)==4
	setappdata(handles.DisplayAxis,'PointerCallback','DisplayPixelCurve');
else
	setappdata(handles.DisplayAxis,'PointerCallback','');
end
setappdata(handles.volViewerFigure,'WaitForClose',false);
set([handles.AcquireBoth handles.AcquireVolume, handles.AcquireSlices handles.AcquireTimeFrames],'enable','off');
i=4;
pos = [];
while i<length(p)
	switch lower(p{i})
		case 'figurename'
			set(handles.volViewerFigure,'Name',p{i+1});
		case 'figurecolor'
			set(handles.volViewerFigure,'Color',p{i+1});
		case 'blockslicecontrol'
			setappdata(handles.DisplayAxis,'BlockSliceControl',logical(p{i+1}));
		case 'framepaneltitle'
			if ischar(p{i+1})
				set(handles.TimeFramePanel,'Title',p{i+1});
			else
				msgID = sprintf('FADSTool:%s:wrongDataType',mfilename);
				error(msgID, 'FramePanelTitle property must be of type character.');
			end
		case 'axisnames'
			titles = p{i+1};
			if length(titles)~=3
				error('Axis names must be a cell of strings with three elements');
			end
			xlabel(handles.DisplayAxis,titles{1});
			ylabel(handles.DisplayAxis,titles{2});
			zlabel(handles.DisplayAxis,titles{3});
			set(handles.XAxisLabel,'string',titles{1});
			set(handles.YAxisLabel,'string',titles{2});
			set(handles.ZAxisLabel,'string',titles{3});
		case 'colorbar'
			if ischar(p{i+1})
				colorbarStatus = p{i+1};
			elseif ~p{i+1}
				colorbarStatus = 'off';
			end
		case 'contours'
			contours = p{i+1};
		case 'extrasurface'
			setappdata(handles.DisplayAxis,'ExtraSurface',p{i+1});
			for ci=1:length(p{i+1})
				contour = p{i+1}{ci};
				if ~isfield(contour,'Type') || isempty(contour.Type) || strcmpi(contour.Type,'mesh')
					h = mesh(handles.DisplayAxis,contour.Y,contour.X,contour.Z);
				else
					h = plot3(contour.Y(:),contour.X(:),contour.Z(:),'parent',handles.DisplayAxis,'Marker','.','LineStyle','-');
				end
				if ischar(contour.LineStyle)
					valid = 'bgrcmykw';
					ii = 1;
					while ii<=length(contour.LineStyle)
						if any(valid==contour.LineStyle(ii))
							ii = ii+1;
						else
							contour.LineStyle = contour.LineStyle(2:end);
						end
					end
				end
				switch get(h,'type')
					case 'surface'
						set(h,'tag',['ExtraSurface' num2str(ci)],'facealpha',0,'edgecolor',contour.LineStyle,'edgelighting','gouraud');
					case 'line'
						set(h,'tag',['ExtraSurface' num2str(ci)],'markerEdgeColor',contour.LineStyle,'markerFaceColor',contour.LineStyle);
				end
				daspect(handles.DisplayAxis,1./getappdata(handles.DisplayAxis,'PixelDimensions'));
			end
		case 'colormap'
			if ischar(p{i+1})
				if lower(p{i+1}(1))=='i'
					colormapInv = true;
					colormapName = p{i+1}(2:end);
				else
					colormapInv = false;
					colormapName = p{i+1};
				end
				set(handles.ColormapInverse,'value',colormapInv);
				cmap=feval(lower(colormapName),128);
				if get(handles.ColormapInverse,'value')
					cmap = flipud(cmap);
				end
				set(handles.volViewerFigure,'colormap',cmap);
			else
				set(handles.volViewerFigure,'Colormap',p{i+1});
				colormapName = 'Custom';
			end
			colorMaps = get(handles.Colormap,'string');
			vi = find(strcmpi(colorMaps,colormapName));
			if isempty(vi)
				vi = 1;
			end
			set(handles.Colormap,'Value',vi);
		case 'time'
			if length(p{i+1}) > 1 && length(s)==4 && length(p{i+1})== s(4)
				setappdata(handles.volViewerFigure,'TimePoints',p{i+1})
			end
			set(handles.TimeOperationPopupmenu,'string',...
				{'Sum';'Average';'Weighted Average';'Integrate';'Max'});
		case 'timeop'
			setPullDownValue(handles.TimeOperationPopupmenu,p{i+1});
		case 'timeunits'
			setappdata(handles.volViewerFigure,'TimeUnits',p{i+1})
		case 'units'
			setappdata(handles.volViewerFigure,'Units',p{i+1})
		case 'pointercallback'
			setappdata(handles.DisplayAxis,'PointerCallback',p{i+1});
		case 'waitforclose'
			setappdata(handles.volViewerFigure,'WaitForClose',logical(p{i+1}));
		case 'position'
			pos = p{i+1};
		case 'pixeldimensions'
			setappdata(handles.DisplayAxis,'PixelDimensions',p{i+1});
			daspect(handles.DisplayAxis,1./p{i+1});
		case 'selectionmode'
			setappdata(handles.volViewerFigure,'SelectionMode',lower(p{i+1}));
			if strcmpi(p{i+1},'Volume')
				set([handles.AcquireSlices, handles.AcquireBoth],'enable','off');
			end
		case 'nargout'
			if p{i+1} > 0 % any output parameters
				selMode = getappdata(handles.volViewerFigure,'SelectionMode');
				if isempty(selMode)
					selMode = 'Auto';
				end
				switch selMode
					case 'volume'
						set(handles.AcquireVolume,'enable','on');
						set([handles.AcquireBoth, handles.AcquireSlices],'enable','off');
					otherwise
						set([handles.AcquireBoth, handles.AcquireVolume, handles.AcquireSlices],'enable','on');
				end
				%need to stop execution until volViewer is closed - since nargout is the last parameter this will override any WaitForClose parameter.
				setappdata(handles.volViewerFigure,'WaitForClose',true);
			end
			if p{i+1} >1
				set(handles.AcquireTimeFrames,'enable','on');
			end
		otherwise
			msgID = sprintf('FADSTool:%s:unknownProperty',mfilename);
			error(msgID, 'There is no ''%s'' property in volViewer', p{i});
	end
	i=i+2;
end
if i==length(p)
	msgID = sprintf('FADSTool:%s:noPropertyValue',mfilename);
	error(msgID, 'No value passed to property ''%s''.', p{i})
end

%% Position
% ========
if ~isempty(pos)
	if ischar(pos)
		movegui(handles.volViewerFigure ,pos);
	else
		set(handles.volViewerFigure,'units','normalized','position',pos);
	end
end
setappdata(handles.volViewerFigure,'DefPos',pos);

%% Colorbar
% ========
updateVolViewerColorBar(handles, colorbarStatus);
for i=1:length(contours)
	arrow(handles.Colorbar,'up',contours(i));
end
setappdata(handles.Colorbar,'ContourValues',[getappdata(handles.Colorbar,'ContourValues') contours]);
FrameSelectionChange(handles.volViewerFigure);
updateContour(handles)

%% Update the main display
% =======================
waitbar(0.5,wbh,'Computing data display. Please wait.');
FrameSelectionChange(handles.volViewerFigure);
waitbar(0.9,wbh,'Applying settings. Please wait.');

if getappdata(handles.volViewerFigure,'WaitForClose')
	set(handles.volViewerFigure,'CloseRequestFcn','set(gcf,''Visible'',''off'');');
else
	set(handles.volViewerFigure,'CloseRequestFcn','try, delete(getappdata(gcf,''ChildObjects'')); catch, end, delete(gcf)');
end
set(handles.volViewerFigure,'WindowButtonMotionFcn',@volViewerPointerMotion,...
							'WindowButtonUpFcn','objectMove');

% Tooltip textbox for labeling value of pointer position
tip = text('Parent',handles.XAxis,'visible','off','fontname','helvetica','fontsize',10,'color','black',...
    'backgroundcolor',[1 1 .8],'edgecolor',[.5 .5 .5],'margin',5,'units','data','tag','Tip');
setappdata(handles.DisplayAxis,'Tip',tip);

% Update the Button Down Callback for all the slices that have already been
% created
h=findobj(get(handles.DisplayAxis,'Children'),'tag','Slice');
for i=1:length(h)
	set([h(i) getappdata(h(i),'Contour')],'buttonDownFcn',[getappdata(handles.DisplayAxis,'PointerCallback') '(gco)']);
end

set(handles.volViewerFigure,'CurrentAxes', handles.DisplayAxis);

close(wbh); % Close the waitbar

set(handles.volViewerFigure,'Pointer','arrow','Visible','On');
