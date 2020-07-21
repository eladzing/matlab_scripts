% Slice Management.  Uses specialized slicomatic slices, not slices
% created with the SLICE command.

function s=updateSlice(handles, X, Y, Z, h)

s=[];

xlim = [getappdata(handles.XAxis,'MinRange') getappdata(handles.XAxis,'MaxRange')];
ylim = [getappdata(handles.YAxis,'MinRange') getappdata(handles.YAxis,'MaxRange')];
zlim = [getappdata(handles.ZAxis,'MinRange') getappdata(handles.ZAxis,'MaxRange')];

if ~isempty(X)
	xi=round(X);
	if xi>=xlim(1) && xi<=xlim(2)
		y = ylim(1) : ylim(2);
		z = zlim(1) : zlim(2);
		data = getappdata(handles.DisplayAxis,'SummedData');
		cdata = reshape(data(xi,y,z),ylim(end)-ylim(1)+1, zlim(end)-zlim(1)+1);
		[xdata, ydata, zdata] = meshgrid(xi, y, z);
		st = 'X';
	else
		return
	end

elseif ~isempty(Y)
	yi=round(Y);
	if yi>=ylim(1) && yi<=ylim(2)
		x = xlim(1) : xlim(2);
		z = zlim(1) : zlim(2);
		data = getappdata(handles.DisplayAxis,'SummedData');
		cdata = reshape(data(x,yi,z),xlim(end)-xlim(1)+1, zlim(end)-zlim(1)+1);
		[xdata, ydata, zdata]=meshgrid(x, yi, z);
		st = 'Y';
	else
		return
	end

elseif ~isempty(Z)
	zi=round(Z);
	if zi>=zlim(1) && zi<=zlim(2)
		x = xlim(1) : xlim(2);
		y = ylim(1) : ylim(2);
		data = getappdata(handles.DisplayAxis,'SummedData');
		cdata = reshape(data(x,y,zi), xlim(end)-xlim(1)+1, ylim(end)-ylim(1)+1)';
		[xdata, ydata, zdata]=meshgrid(x, y, zi);
		st = 'Z';
	else
		return
	end
else
	error('No slice information was passed to updateSlice.');
end

cdata=squeeze(cdata);
xdata=squeeze(xdata);
ydata=squeeze(ydata);
zdata=squeeze(zdata);

if ishandle(handles.Colorbar)
	levels = getappdata(handles.Colorbar,'ContourValues');
	if ~isempty(levels)
		if length(levels)==1
			levels = [levels levels];
		end
		switch st
			case 'X'
				c = contours(zdata,ydata,cdata,levels);
			case 'Y'
				c = contours(zdata,xdata,cdata,levels);
			case 'Z'
				c = contours(xdata,ydata,cdata,levels);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		newvertices = [];
		newfaces = {};
		longest = 1;
		ccdata = [];

		i = 1;
		while(i < size(c,2))
			z_level = c(1,i);
			npoints = c(2,i);
			longest=max(longest,npoints);

			cxdata = c(1,i+1:i+npoints);
			cydata = c(2,i+1:i+npoints);

			switch st
				case 'X'
					lzdata = xi*ones(1,npoints);
					vertices = [lzdata', cydata', cxdata'];
				case 'Y'
					lzdata = yi*ones(1,npoints);
					vertices = [cydata', lzdata', cxdata'];
				case 'Z'
					lzdata = zi*ones(1,npoints);
					vertices = [cxdata', cydata', lzdata'];
			end

			faces = size(newvertices,1)+(1:npoints);

			newvertices = [ newvertices ; vertices ];
			newfaces{end+1} = faces;

			tcdata =  z_level*ones(npoints,1);

			ccdata = [ ccdata; tcdata ]; % need to be same size as faces

			i = i+npoints+1;
		end

		% Glom a NAN on the end for loop-breaking
		newvertices = [ newvertices ; nan nan nan ];
		ccdata = [ ccdata ; nan ];

		vertmax = size(newvertices,1);

		% Fix up FACES, which is a cell array.
		faces = [];
		for i = 1:size(newfaces,2)
			faces = [ faces; newfaces{i} ones(1,longest-size(newfaces{i},2))*vertmax vertmax ];
		end
	else % NO CONTOURS
		newvertices=[];
		faces=[];
		ccdata=[];
	end
else % NO COLORBAR OBJECT FOUND
	newvertices=[];
	faces=[];
	ccdata=[];
end

if nargin == 5 && ~isempty(h) % Recycle the old slice
	set(h,'cdata',cdata,'alphadata',cdata, 'xdata',xdata, ...
		'ydata',ydata, 'zdata',zdata);
	s=h;
	if strcmp(get(s,'facec'),'texturemap')
		textureizeslice(s,'on');
	end
	set(getappdata(h,'Contour'),'vertices',newvertices,'faces',faces,'facevertexcdata',ccdata); % update the contour patch

else % Create a new slice
	h=getappdata(handles.volViewerFigure,'SliceContextMenu');
	s=surface('Parent',handles.DisplayAxis,'tag','Slice',...
		'cdata',cdata,'xdata',xdata, 'ydata',ydata, 'zdata',zdata,...
		'alphadata',abs(cdata),'alphadatamapping','scaled',...
		'facelighting','none',...
		'uicontextmenu',h.sliceMenu,...
		'buttonDownFcn',[getappdata(handles.DisplayAxis,'PointerCallback') '(gco)']);
	setappdata(s,'SliceType',st);
	% Use default color and alpha setting
	color = getappdata(handles.volViewerFigure,'DefaultColor');
	if strcmp(color,'faceted')
		set(s,'facecolor','flat','edgec','k');
	else %flat, interp, texture
		set(s,'facecolor',color,'edgec','n');
	end
	alpha = getappdata(handles.volViewerFigure,'DefaultAlpha');
	if strcmp(alpha,'none')
		set(s,'facealpha',1);
	else   %flat, interp, texture
		set(s,'facealpha',alpha);
	end
	ch = patch('Tag','Contour',...
		'Parent',handles.DisplayAxis,...
		'vertices',newvertices,'faces',faces,'facevertexcdata',ccdata,...
		'uicontextmenu',h.sliceMenu,...
		'buttonDownFcn',getappdata(handles.DisplayAxis,'PointerCallback')); % create a patch for the contour to go with the slice
	setappdata(ch,'Slice',s); % tie the slice to the contour
	setappdata(s,'Contour',ch);
end