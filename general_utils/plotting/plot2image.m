function [ h_axis ] = plot2image( axis_obj, resolution)
% PLOT2IMAGE Generates an axis object copied from the input axis, but
% with the plot contents reproduced as an image. The idea is to eb able to
% save figure with shaded patch objects as vector images. 
%
% Syntax:
% PLOT2IMAGE(axis_object) will create a new figure with an axis copied
% from the axis referenced by the axis_object handle, but with the content
% of the plot converted to an image. 
% PLOT2IMAGE(axis_object, resolution) defines the resolution of the
% generated image (used in the "print" command) in dots per inch (dpi). 
% Default value is 100. 
% h = PLOT2IMAGE(...) returns the handle of the new axis object. 
% 
% Written by Edden Gerber, lab of Leon Y. Deouell, Hebrew University of
% Jerusalem, 2015. 
%

% Optional input arguments
if nargin < 2
    resolution = 100;
end

% Create new figure
h_fig = figure;

% Get tick positions
Xlim = axis_obj.XLim;
Ylim = axis_obj.YLim;
xtick = axis_obj.XTick;
ytick = axis_obj.YTick;

% Get labels
xlab = axis_obj.XLabel;
ylab = axis_obj.YLabel;
ti = axis_obj.Title;

% Make axis fill figure so only it will be captured by the print command
orig_pos = axis_obj.Position;
axis_obj.Position = [0 0 1 1];
drawnow;

% get rid of tick marks and borders before taking snapshot
axis off

% get snapshot and flip y axis (to make origin bottom left corner)
img = print(axis_obj.Parent,'-RGBImage',['-r' num2str(resolution)]);
img = img(end:-1:1,:,:);

% restore original axis
axis_obj.Position = orig_pos;
axis on

% Get image dimensions
img_dim = [size(img,1) size(img,2)];

% Make new axes
new_x_axis = linspace(Xlim(1),Xlim(2),img_dim(2));
new_y_axis = linspace(Ylim(1),Ylim(2),img_dim(1));

% Display image
image(new_x_axis,new_y_axis, img); axis xy;
h_axis = gca;
h_axis.XTick = xtick;
h_axis.YTick = ytick;
h_axis.XLabel = copy(xlab); % if the label object is not copied, it will be taken from the orignial... 
h_axis.YLabel = copy(ylab);
h_axis.Title = copy(ti);

end

