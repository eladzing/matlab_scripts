% Set up volViwer menus. This includes the general GUI menus and context
% menus specific to certain objects. The menu handles are recorded as
% appdata in the DisplayAxis object.

% By Ran Klein, University of Ottawa Heart Institute, 2-Nov-2005

function volViewerMenus(p)


if ~ischar(p)
	% Main Figure Menu
	handles = p;
	fig = handles.volViewerFigure;
	set(fig,'menubar','none'); %

	if ~isempty(findobj(get(fig,'children'),'flat','tag','FileMenu'))
		return;
	end
	
	% File menu
	h.fileMenu = uimenu(fig,'label','File','Tag','FileMenu');
	h.figCopy = uimenu(h.fileMenu, 'label', 'Copy figure','callback', 'volViewerMenus Copy');
	h.figPrint  = uimenu(h.fileMenu,'label','Print...','callback','volViewerMenus Print');
	h.fsaveprefs = uimenu(h.fileMenu,'label','Save preferences','callback','volViewerMenus SavePrefs');
	h.fexit = uimenu(h.fileMenu, 'label', 'Close','callback','set(gcf,''Visible'',''off'');',...
		'separator','on');

	% Controls Menu
	h.controlMenu = uimenu(fig,'label','Controls', 'callback',@controlmenu,'Tag','ControlMenu');
	h.camtoolbar = uimenu(h.controlMenu,'label','Camera toolbar','callback', 'volViewerMenus CameraToolbar');
	h.dclabels= uimenu(h.controlMenu','label','Tick Labels','callback','volViewerMenus ControlLabels','Separator','On');
	h.dcvis  = uimenu(h.controlMenu','label','Control Bars','callback','volViewerMenus ControlVisible');
	h.colorbar = uimenu(h.controlMenu','label','Color Bar','callback','volViewerMenus Colorbar');
	h.coordmode  = uimenu(h.controlMenu','label','Coordinate Mode','callback',@coordmodemenu,'Tag','CoordModeMenu');
	h.xyz  = uimenu(h.coordmode','label','Matrix (XYZ)','callback','volViewerMenus XYZ');
	h.sct  = uimenu(h.coordmode','label','Subject (SCT)','callback','volViewerMenus SCT');
	h.hvsh  = uimenu(h.coordmode','label','Camera (HVSh)','callback','volViewerMenus HVSh');

	% Default for new slices menu
	h.defMenu = uimenu(fig,'label','Slice Defaults', 'callback', @defaultmenu,'Tag','DefaultMenu');
	h.dfacet  = uimenu(h.defMenu,'label','Slice Color Faceted','callback','volViewerMenus DefaultFaceted');
	h.dflat   = uimenu(h.defMenu,'label','Slice Color Flat',   'callback','volViewerMenus DefaultFlat');
	h.dinterp = uimenu(h.defMenu,'label','Slice Color Interp', 'callback','volViewerMenus DefaultInterp');
	h.dtex    = uimenu(h.defMenu,'label','Slice Color Texture','callback','volViewerMenus DefaultTexture');
	h.dcnone  = uimenu(h.defMenu,'label','Slice Color None','callback','volViewerMenus DefaultColorNone');
	h.dtnone  = uimenu(h.defMenu,'label','Slice Transparency None','callback','volViewerMenus DefaultTransNone','separator','on');
	h.dtflat  = uimenu(h.defMenu,'label','Slice Transparency Flat','callback','volViewerMenus DefaultTransFlat');
	h.dtinterp= uimenu(h.defMenu,'label','Slice Transparency Interp','callback','volViewerMenus DefaultTransInterp');
	h.dttex   = uimenu(h.defMenu,'label','Slice Transparency Texture','callback','volViewerMenus DefaultTransTexture');

	% Set default values
	setappdata(handles.volViewerFigure,'DefaultColor','texture');
	setappdata(handles.volViewerFigure,'DefaultAlpha','texture');
	setappdata(handles.volViewerFigure,'MaskCoordMode',{'X' 'Y' 'Z'});
	LoadUserPreferences(handles);
	
	% Set props for all slices menu
	h.allMenu = uimenu(fig,'label','All Slices','Tag','AllSlicesMenu');
	uimenu(h.allMenu,'label','Color Faceted','callback','volViewerMenus AllFacet');
	uimenu(h.allMenu,'label','Color Flat','callback','volViewerMenus AllFlat');
	uimenu(h.allMenu,'label','Color Interp','callback','volViewerMenus AllInterp');
	uimenu(h.allMenu,'label','Color Texture','callback','volViewerMenus AllTex');
	uimenu(h.allMenu,'label','Color None','callback','volViewerMenus AllNone');
	uimenu(h.allMenu,'label','Transparency None','callback','volViewerMenus AllTNone','separator','on');
	uimenu(h.allMenu,'label','Transparency .5','callback','volViewerMenus AllTP5');
	uimenu(h.allMenu,'label','Transparency Flat','callback','volViewerMenus AllTFlat');
	uimenu(h.allMenu,'label','Transparency Interp','callback','volViewerMenus AllTInterp');
	uimenu(h.allMenu,'label','Transparency Texture','callback','volViewerMenus AllTTex');

	% Setup Help style options
	h.helpmenu = uimenu(fig,'label','Help','Tag','HelpMenu');
	uimenu(h.helpmenu,'label','Help','callback','system ''volViewer Manual.pdf''');
% 	uimenu(h.helpmenu,'label','Check for Updates','callback','web http://www.site.uottawa.ca/~rklein/Personal/volumeViewer.htm');
	uimenu(h.helpmenu,'label','About Author','callback','web http://www.ottawaheart.ca/misc/ran-klein.htm');

	setappdata(handles.volViewerFigure,'Menus',h);


	% Context Menus
	% =============
	% Slice Context Menu
	h=[];
	h.sliceMenu=uicontextmenu('Parent',fig,'callback', @slicecontextmenu,'Tag','SliceMenu');
	h.vistog = uimenu(h.sliceMenu,'label','Visible','callback','volViewerMenus ToggleVisible');
	h.sliceDelete = uimenu(h.sliceMenu,'label','Delete','callback','volViewerMenus DeleteSlice');
	h.smcolorm  = uimenu(h.sliceMenu,'label','Color','separator','on');
	h.smfacet   = uimenu(h.smcolorm,'label','Color Faceted','callback','volViewerMenus Setfaceted');
	h.smflat    = uimenu(h.smcolorm,'label','Color Flat','callback','volViewerMenus SetFlat');
	h.sminterp  = uimenu(h.smcolorm,'label','Color Interp','callback','volViewerMenus SetInterp');
	h.smtex     = uimenu(h.smcolorm,'label','Color Texture','callback','volViewerMenus SetTexture');
	h.smnone    = uimenu(h.smcolorm,'label','Color None','callback','volViewerMenus SetNone');
	h.smtransm  = uimenu(h.sliceMenu,'label','Transparency');
	h.smtnone   = uimenu(h.smtransm,'label','Transparency None','callback','volViewerMenus SetAlphaNone');
	h.smtp5     = uimenu(h.smtransm,'label','Transparency .5','callback','volViewerMenus SetAlphaPoint5');
	h.smtflat   = uimenu(h.smtransm,'label','Transparency Flat','callback','volViewerMenus SetAlphaFlat');
	h.smtinterp = uimenu(h.smtransm,'label','Transparency Interp','callback','volViewerMenus SetAlphaInterp');
	h.smttex    = uimenu(h.smtransm,'label','Transparency Texture','callback','volViewerMenus SetAlphaTexture');
	h.showsurf = uimenu(h.sliceMenu,'label','Show Surface','Separator','on','callback','volViewerMenus ShowSurface');

	setappdata(handles.volViewerFigure,'SliceContextMenu',h);

else % called by a callback function
	handles = guidata(gcbo);
	h = getappdata(handles.volViewerFigure,'Menus');
	switch p
		case 'Copy'
			fig = figure;
			h = copyobj(handles.DisplayAxis,fig, 'legacy');
			set(h,'pos',[.1 .1 .8 .8],'units','normalized');
			if ~isempty(handles.Colorbar)
				colorbar('peer',h);
			end
			print(fig,'-dmeta');
		case 'Print'
			newf=figure('visible','off','renderer',get(handles.volViewerFigure,'renderer'));
			h = copyobj(handles.DisplayAxis, newf, 'legacy');
			set(gca,'pos',[.1 .1 .8 .8],'units','normalized')
			if ~isempty(handles.Colorbar)
				colorbar('peer',h);
			end
			printdlg(newf);
			close(newf);
		case 'SavePrefs'
			%extract only preferences that need to be sticky
			prefs.defcolor = getappdata(handles.volViewerFigure,'DefaultColor');
			prefs.defalpha = getappdata(handles.volViewerFigure,'DefaultAlpha');
			prefs.maskcoordmode = getappdata(handles.volViewerFigure,'MaskCoordMode');
			prefs.camtoolbar_checked = cameratoolbar('getvisible');
			prefs.ticklabels = get(handles.XAxis,'xticklabelmode');
			prefs.colormap = get(handles.Colormap,'string');
			prefs.colormap = prefs.colormap{get(handles.Colormap,'value')};
			if strcmpi(prefs.colormap,'custom')
				prefs.colormap = get(handles.volViewerFigure,'Colormap');
			end
			prefs.colormapInv = get(handles.ColormapInverse,'value');
			fileName = UserStickyPrefsFileName;
			save(fileName,'prefs')

		% Control
		% -------
		case 'CameraToolbar'
			cameratoolbar('Toggle');
		case 'ControlLabels'
			l = get(handles.XAxis,'xticklabel');
			if isempty(l)
				set([handles.XAxis],'xticklabelmode','auto');
				set([handles.YAxis, handles.ZAxis],'yticklabelmode','auto');
			else
				set([handles.XAxis],'xticklabel',[]);
				set([handles.YAxis, handles.ZAxis],'yticklabel',[]);
			end
		case 'ControlVisible'
			objects = findobj([handles.XAxis, handles.YAxis, handles.ZAxis,...
				handles.XAxisLabel, handles.YAxisLabel, handles.ZAxisLabel]);
			if strcmp(get(handles.XAxis,'Visible'),'on')
				set(objects,'visible','off');
			else
				set(objects,'visible','on');
			end
		case 'Colorbar'
			if strcmpi(get(gcbo,'Checked'),'On')
				updateVolViewerColorBar(handles,'Off');
				set(gcbo,'Checked','Off');
			else
				updateVolViewerColorBar(handles,'On');
				set(gcbo,'Checked','On');
			end
		case 'XYZ', setappdata(handles.volViewerFigure,'MaskCoordMode',{'X' 'Y' 'Z'});
		case 'SCT', setappdata(handles.volViewerFigure,'MaskCoordMode',{'S' 'C' 'T'});
		case 'HVSh', setappdata(handles.volViewerFigure,'MaskCoordMode',{'H' 'V' 'Sh'});
		% Defaults
		% --------
		case 'DefaultFaceted', setappdata(handles.volViewerFigure,'DefaultColor','faceted');
		case 'DefaultFlat', setappdata(handles.volViewerFigure,'DefaultColor','flat');
		case 'DefaultInterp', setappdata(handles.volViewerFigure,'DefaultColor','interp');
		case 'DefaultTexture', setappdata(handles.volViewerFigure,'DefaultColor','texture'); 
			if strcmp(getappdata(handles.volViewerFigure,'DefaultAlpha'),'flat') || strcmp(getappdata(handles.volViewerFigure,'DefaultAlpha'),'interp')
				setappdata(handles.volViewerFigure,'DefaultAlpha','texture');
			end
		case 'DefaultColorNone', setappdata(handles.volViewerFigure,'DefaultColor','none');
		case 'DefaultTransNone', setappdata(handles.volViewerFigure,'DefaultAlpha','none');
		case 'DefaultTransFlat', setappdata(handles.volViewerFigure,'DefaultAlpha','flat');
		case 'DefaultTransInterp', setappdata(handles.volViewerFigure,'DefaultAlpha','interp');
		case 'DefaultTransTexture', setappdata(handles.volViewerFigure,'DefaultAlpha','texture'); 
			setappdata(handles.volViewerFigure,'DefaultColor','texture');


		% All Slices
		% ----------
		case 'AllFacet'
			s=findobj(handles.DisplayAxis,'type','surface','tag','Slice');
			set(s,'facec','flat','edgec','k');
			textureizeslice(s,'off');
		case 'AllFlat'
			s=findobj(handles.DisplayAxis,'type','surface','tag','Slice');
			set(s,'facec','flat','edgec','none');
			textureizeslice(s,'off');
		case 'AllInterp'
			s=findobj(handles.DisplayAxis,'type','surface','tag','Slice');
			set(s,'facec','interp','edgec','none');
			textureizeslice(s,'off');
		case 'AllTex'
			s=findobj(handles.DisplayAxis,'type','surface','tag','Slice');
			set(s,'facec','texturemap','edgec','none');
			textureizeslice(s,'on');
		case 'AllNone'
			s=findobj(handles.DisplayAxis,'type','surface','tag','Slice');
			set(s,'facec','none','edgec','none');
			textureizeslice(s,'off');
		case 'AllTNone'
			s=findobj(handles.DisplayAxis,'type','surface','tag','Slice');
			set(s,'facea',1);
			textureizeslice(s,'off');
		case 'AllTP5'
			s=findobj(handles.DisplayAxis,'type','surface','tag','Slice');
			set(s,'facea',.5);
			textureizeslice(s,'off');
		case 'AllTFlat'
			s=findobj(handles.DisplayAxis,'type','surface','tag','Slice');
			set(s,'facea','flat');
			textureizeslice(s,'off');
		case 'AllTInterp'
			s=findobj(handles.DisplayAxis,'type','surface','tag','Slice');
			set(s,'facea','interp');
			textureizeslice(s,'off');
		case 'AllTTex'
			s=findobj(handles.DisplayAxis,'type','surface','tag','Slice');
			set(s,'facea','texturemap');
			textureizeslice(s,'on');
		otherwise
		% Context Menus
		% -------------
		ah = gco;
		if strcmp(get(ah,'Tag'),'SliceArrow')
			sh = getappdata(ah,'Slice');
		elseif strcmp(get(ah,'Tag'),'Slice')
			sh = ah;
			ah = getappdata(sh,'Arrow');
		elseif strcmp(get(ah,'Tag'),'Surface')
			sh = ah;
			ah = getappdata(sh,'Arrow');
		else
			error('Unexpected object called sliceContextMenu or unsupported menu operation ''%s'' detected by %s.', p, mfilename);
		end
		switch p
			case 'ToggleVisible'
				if strcmp(get(sh,'Tag'),'Slice')
					sh = [sh getappdata(sh,'Contour')];
					if strcmp(get(sh,'visible'),'on')
						set(sh,'visible','off');
						set(ah,'facecolor','y');
					else
						set(sh,'visible','on');
						set(ah,'facecolor','g');
					end
				else
					if isequal(get(ah,'FaceColor'),[0 1 0]) % only if visible
						set(ah,'facecolor','y');
					else
						set(ah,'facecolor','g');
					end
					updateContourVals(handles.Colorbar)
					updateContour(handles)
				end
			case 'DeleteSlice'
				if get(ah,'Parent')==handles.Colorbar
					if ~isempty(getappdata(ah,'Surface'))
						delete(getappdata(ah,'Surface'));
					end
					delete(ah);
					values = [];
					h = findobj(get(handles.Colorbar,'Children'),'tag','SliceArrow');
					for i=1:length(h)
						values = [values getappdata(h(i),'CenterPosition')];
					end
					setappdata(handles.Colorbar,'ContourValues',values);
					updateContour(handles);
				else
					delete(getappdata(sh,'Contour'));
					delete(ah); delete(sh);
					h=findobj([handles.XAxis, handles.YAxis, handles.ZAxis],'Type','patch','Tag','SliceArrow');
					if isempty(h) % Any other slice left?
						set(handles.AcquireSlices,'Enable','Off');
					end
				end
			case 'Setfaceted'
				set(sh,'edgec','k','facec','flat');
				if ischar(get(sh,'facea')) && strcmp(get(sh,'facea'),'texturemap')
					set(sh,'facea','flat');
				end
				textureizeslice(sh,'off');
			case 'SetFlat'
				set(sh,'edgec','n','facec','flat');
				if ischar(get(sh,'facea')) && strcmp(get(sh,'facea'),'texturemap')
					set(sh,'facea','flat');
				end
				textureizeslice(sh,'off');
			case 'SetInterp'
				set(sh,'edgec','n','facec','interp');
				if ischar(get(sh,'facea')) && strcmp(get(sh,'facea'),'texturemap')
					set(sh,'facea','interp');
				end
				textureizeslice(sh,'off');
			case 'SetTexture'
				set(sh,'facecolor','texture','edgec','none');
				if ischar(get(sh,'facea'))
					set(sh,'facealpha','texturemap');
				end
				textureizeslice(sh,'on');
			case 'SetNone'
				set(sh,'facecolor','none','edgec','none');
				textureizeslice(sh,'off');
			case 'SetAlphaNone'
				if get(ah,'Parent')==handles.Colorbar
					set(getappdata(ah,'Surface'),'FaceAlpha',1);
				else
					set(sh,'facealpha',1);
				end
			case 'SetAlphaPoint5'
				if get(ah,'Parent')==handles.Colorbar
					set(getappdata(ah,'Surface'),'FaceAlpha',0.5);
				else
					set(sh,'facealpha',.5);
				end
			case 'SetAlphaFlat'
				set(sh,'facealpha','flat');
				if ischar(get(sh,'facec')) && strcmp(get(sh,'facec'),'texturemap')
					set(sh,'facecolor','flat');
					textureizeslice(sh,'off');
				end
			case 'SetAlphaInterp'
				set(sh,'facealpha','interp');
				if ischar(get(sh,'facec')) && strcmp(get(sh,'facec'),'texturemap')
					set(sh,'facecolor','interp');
					textureizeslice(sh,'off');
				end
			case 'SetAlphaTexture'
				set(sh,'facealpha','texturemap');
				if ischar(get(sh,'facec'))
					set(sh,'facecolor','texturemap');
					textureizeslice(sh,'on');
				end
			case 'ShowSurface'
				if isempty(getappdata(ah,'Surface'))
					updateSurface(ah,handles);
				else
					delete(getappdata(ah,'Surface'));
					setappdata(ah,'Surface',[]);
				end
			otherwise
				error('Unsupported menu operation ''%s'' detected by %s.', p, mfilename);
		end
	end
end

% Function that updates the setting in the control menu when it is opened
function controlmenu(hObject, eventdata, handles)

handles=guidata(hObject);
h=getappdata(handles.volViewerFigure,'Menus');

if cameratoolbar(handles.volViewerFigure,'getvisible')
	set(h.camtoolbar,'checked','on');
else
	set(h.camtoolbar,'checked','off');
end

if isempty(get(handles.XAxis,'xticklabel'))
	set(h.dclabels,'checked','off');
else
	set(h.dclabels,'checked','on');
end

set(h.dcvis,'checked',get(handles.XAxis,'visible'));
set(h.colorbar,'checked',get(handles.Colorbar,'Visible'));

% Function that updates the setting in the coordinate mode menu when it is opened
function coordmodemenu(hObject, eventdata, handles)

handles=guidata(hObject);
h=getappdata(handles.volViewerFigure,'Menus');

set([h.xyz h.sct h.hvsh],'checked','off');
m=getappdata(handles.volViewerFigure,'MaskCoordMode');
switch m{1}
	case 'X'
		set(h.xyz,'checked','on');
	case 'S'
		set(h.sct,'checked','on');
	case 'H'
		set(h.hvsh,'checked','on');
	otherwise
		error('Unidentified mask coordinate mode (%s,%s,%s) detected in coordmodemenu.', m{1}, m{2}, m{3});
end


% Function that updates the display in the default settings menu when it is
% open
function defaultmenu(fig, action)

handles=guidata(gcbo);
h=getappdata(handles.volViewerFigure,'Menus');

set([h.dfacet h.dflat h.dinterp h.dtex h.dtnone h.dtflat h.dtinterp h.dttex],...
	'checked','off');
switch getappdata(handles.volViewerFigure,'DefaultColor')
	case 'faceted', set(h.dfacet,'checked','on');
	case 'flat', set(h.dflat,'checked','on');
	case 'interp', set(h.dinterp,'checked','on');
	case 'texture', set(h.dtex,'checked','on');
	case 'none', set(h.dcnone,'checked','on');
end
switch getappdata(handles.volViewerFigure,'DefaultAlpha')
	case 'none', set(h.dtnone,'checked','on');
	case 'flat', set(h.dtflat,'checked','on');
	case 'interp', set(h.dtinterp,'checked','on');
	case 'texture', set(h.dttex,'checked','on');
end

%% Context menu state for slices
function slicecontextmenu(fig,action)
ah = gco;
if strcmp(get(ah,'Tag'),'SliceArrow')
	sh = getappdata(ah,'Slice');
elseif strcmp(get(ah,'Tag'),'Slice')
	sh = ah;
	ah = getappdata(sh,'Arrow');
elseif strcmp(get(ah,'Tag'),'Surface')
	sh = ah;
	ah = getappdata(sh,'Arrow');
else
	error('Unexpected object called sliceContextMenu.');
end
	
handles=guidata(ah);
h=getappdata(handles.volViewerFigure,'SliceContextMenu');

if get(ah,'parent')==handles.Colorbar
	% specific items for contour context-menus
	set(h.smcolorm,'visible','off');
	set(h.showsurf,'visible','on');
	
	if isequal(get(ah,'FaceColor'),[0 1 0]) % only if visible
		set(h.vistog,'checked','on');
	else
		set(h.vistog,'checked','off');
	end
	if ~isempty(getappdata(ah,'Surface'))
		set(h.showsurf,'checked','on');
		set(h.smtransm,'visible','on');
		set([h.smtflat, h.smtinterp, h.smttex],'visible','off');
		set([h.smtnone h.smtp5],'checked','off');
		if get(getappdata(ah,'Surface'),'FaceAlpha')==1
			set(h.smtnone,'checked','on');
		else
			set(h.smtp5,'checked','on');
		end
	else
		set(h.showsurf,'checked','off');
		set(h.smtransm,'visible','off');
	end
else
	% specific menu items for slice contextmenus
	set([h.smcolorm, h.smtransm, h.smtflat, h.smtinterp, h.smttex],'visible','on');
	set(h.showsurf,'visible','off');
	% reset all checkmarks - they will be added as needed
	set([h.smfacet h.smflat h.sminterp h.smtex h.smtnone h.smtp5 ...
		h.smtflat h.smtinterp h.smttex h.smnone],'checked','off');

	set(h.vistog,'checked',get(sh,'visible'));

	if isequal(get(sh,'edgec'),[0 0 0])
		set(h.smfacet,'checked','on');
	elseif strcmp(get(sh,'facec'),'flat')
		set(h.smflat,'checked','on');
	end
	if strcmp(get(sh,'facec'),'interp')
		set(h.sminterp,'checked','on');
	end
	if strcmp(get(sh,'facec'),'texturemap')
		set(h.smtex,'checked','on');
	end
	if strcmp(get(sh,'facec'),'none')
		set(h.smnone,'checked','on');
	end
	if isequal(get(sh,'facea'),1)
		set(h.smtnone,'checked','on');
	end
	if isequal(get(sh,'facea'),.5)
		set(h.smtp5,'checked','on');
	end
	if strcmp(get(sh,'facea'),'flat')
		set(h.smtflat,'checked','on');
	end
	if strcmp(get(sh,'facea'),'interp')
		set(h.smtinterp,'checked','on');
	end
	if strcmp(get(sh,'facea'),'texturemap')
		set(h.smttex,'checked','on');
	end
end

% Function that loads the saved defaults
function LoadUserPreferences(handles)

fileName = UserStickyPrefsFileName;

%override particular field values (if file exists)
if exist(fileName,'file')
	try
		load(fileName)
		if prefs.camtoolbar_checked
			cameratoolbar(handles.volViewerFigure,'show');
		else
			cameratoolbar(handles.volViewerFigure,'hide');
		end
		setappdata(handles.volViewerFigure,'DefaultColor',prefs.defcolor);
		setappdata(handles.volViewerFigure,'DefaultAlpha',prefs.defalpha);
		setappdata(handles.volViewerFigure,'MaskCoordMode',prefs.maskcoordmode);

		if strcmp('auto',prefs.ticklabels)
			set([handles.XAxis],'xticklabelmode','auto');
			set([handles.YAxis handles.ZAxis],'yticklabelmode','auto');
		else
			set([handles.XAxis],'xticklabel',[]);
			set([handles.YAxis handles.ZAxis],'yticklabel',[]);
		end
		
		set(handles.ColormapInverse,'value',prefs.colormapInv);
		if ischar(prefs.colormap)
			ncolors = length(get(handles.volViewerFigure,'Colormap'));
			if isempty(ncolors) || ncolors<2
				ncolors = 64;
			end
			cmap=eval([lower(prefs.colormap) '(' num2str(ncolors) ')']);
			if get(handles.ColormapInverse,'value')
				cmap = flipud(cmap);
			end
			set(handles.volViewerFigure,'colormap',cmap);
		else
			set(handles.volViewerFigure,'Colormap',prefs.colormap);
			prefs.colormap = 'Custom';
		end
		colorMaps = get(handles.Colormap,'string');
		set(handles.Colormap,'Value',find(strcmpi(colorMaps,prefs.colormap),1));
		
	catch
		disp('Preferences could not be loaded. Using default values.')
	end
end


%----------------------------------------------------------------------
function fileName = UserStickyPrefsFileName
localPath = fileparts(which(mfilename));
fileName = fullfile(localPath,'volViewer.mat');
