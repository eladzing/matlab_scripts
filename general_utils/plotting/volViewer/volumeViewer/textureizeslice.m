% Convert a regular slice into a texture map slice, or a texture
% slice into a regular slice.


function textureizeslice(slice, onoff)

for k=1:numel(slice)

	d=getappdata(slice(k),'textureoptimizeations');

	switch onoff
		case 'on'
			d.xdata=get(slice(k),'xdata');
			d.ydata=get(slice(k),'ydata');
			d.zdata=get(slice(k),'zdata');
			setappdata(slice(k),'textureoptimizeations',d);
			if max(size(d.xdata)==1)
				nx=[d.xdata(1) d.xdata(end)];
			else
				nx=[d.xdata(1,1)   d.xdata(1,end);
					d.xdata(end,1) d.xdata(end,end)];
			end
			if max(size(d.ydata)==1)
				ny=[d.ydata(1) d.ydata(end)];
			else
				ny=[d.ydata(1,1)   d.ydata(1,end);
					d.ydata(end,1) d.ydata(end,end)];
			end
			if max(size(d.zdata)==1)
				nz=[d.zdata(1) d.zdata(end)];
			else
				nz=[d.zdata(1,1)   d.zdata(1,end);
					d.zdata(end,1) d.zdata(end,end)];
			end
			set(slice(k),'xdata',nx, 'ydata', ny, 'zdata', nz,...
				'facec','texturemap');
			if ischar(get(slice(k),'facea'))
				set(slice(k),'facea','texturemap');
			end
			if ischar(get(slice(k),'facec'))
				set(slice(k),'facec','texturemap');
			end
		case 'off'
			if ~isempty(d)
				set(slice(k),'xdata',d.xdata,'ydata',d.ydata,'zdata',d.zdata);
				setappdata(slice(k),'textureoptimizeations',[]);
			end
			if ischar(get(slice(k),'facea')) && strcmp(get(slice(k),'facea'),'texturemap')
				set(slice(k),'facea','flat');
			end
			if ischar(get(slice(k),'facec')) && strcmp(get(slice(k),'facec'),'texturemap')
				set(slice(k),'facec','flat');
			end
	end
end
