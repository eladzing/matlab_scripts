% volViewerDemo - load demo data and call volViewer.
%
% See comments below about modes of operation in the demo.
% Type 'help volViewer' to get help on the volViewer function usage.
%
% See also: volViewer, DisplayPixelCurve

% By Ran Klein, University of Ottawa Heart Institute, 2014-01-15


%% Change mode of operation:
% false - Do not stop and wait. Just generate display and finish execution.
% true -  Stop and wait for user to close View4D. Returns the time frames
%         selected to generate the last display.
waitForClose = true;

%% Load demo data
load('volViewerDemoData.mat')

%% Call 4DViewer with all options
[ROIs,frames] = volViewer(vol,... % 4D volume data
	'(39,:,:)',...
	'13-15',... % time frames to display at startup
	'PixelDimensions',[xdim, ydim, zdim],... % pixel spacial dimension sizes
	'Time',time,...
	'TimeUnits',timeUnits,... frame times and units
	'TimeOp','Sum',...	
	'Contours',[2.7e8 5e8],...
	'Units',uptakeUnits,... image units name
	'ExtraSurface',ROIdata,... % contours/surface data to add to the slices
	'Position','west',... figure position
	'FigureName','volViewer Demo',... % Name of the figure
	'FigureColor','w',... % Background color of the figure
	'FramePanelTitle','Frames',...	
	'AxisNames',{'Coronal','Sagittal','Transaxial'},... axis labels   
	'Colormap','iHotMetal',... The image colormap
	'BlockSliceControl',false,...
	'WaitForClose',waitForClose,...
	'PointerCallback','DisplayPixelCurve'); % The callback function used to generate a live display when a new intersection pixel is selected

disp('Done View4D call');

%% Display results
if waitForClose
	disp(['ROI ' ROIs ' were selected using volViewer']);
	disp(['Frames ' frames ' was selected using volViewer']);
end