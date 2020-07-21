% updateContour - update display of all the contours.

% By Ran Klein, University of Ottawa Heart Institute, 20-Oct-2005

function updateContour(handles)

h = findobj(get(handles.DisplayAxis,'Children'),'Tag','Slice');
for i=1:length(h)
	st = getappdata(h(i),'SliceType');
	value = getappdata(getappdata(h(i),'Arrow'),'CenterPosition');
	switch st
		case 'X'
			updateSlice(handles,value,[],[],h(i));
		case 'Y'
			updateSlice(handles,[],value,[],h(i));
		case 'Z'
			updateSlice(handles,[],[],value,h(i));
		otherwise
			error('Unrecognized slice type');
	end
end