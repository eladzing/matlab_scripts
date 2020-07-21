% updateContourVals - Updates the values for which contour lines are drawn
% by looking at the location of all the slice arrows on the colorbar. The
% values are stored as an appdata item attached to the colorbar object.
%
% Usage:
% ------
% updateContourVals(ph) where ph is the Colorbar handle.

% By Ran Klein, University of Ottawa Heart Institute, 20-Oct-2005

function updateContourVals(ph)

values = [];
h = findobj(get(ph,'Children'),'tag','SliceArrow');
for i=1:length(h)
	if isequal(get(h(i),'FaceColor'),[0 1 0]) % only if visible
		values = [values getappdata(h(i),'CenterPosition')];
	end
end
setappdata(ph,'ContourValues',values);