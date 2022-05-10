function   mkmapStars(varargin )
%MKMAP plotting the TNG objects
%   Plotting the gas cells from the TNG by mapping to a uniform grid
%% defuals and globals
global DEFAULT_PRINTOUT_DIR
global simDisplayName
global illUnits
global cosmoStruct
units;

projectiontags={'YZ' 'XZ' 'XY'};
bartagOveride='';
printoutdir=DEFAULT_PRINTOUT_DIR;
printtag='badtag';

Ngrid=256;
plotproj=[0 0 0]; % plot all three projections
dilute=4;
streamDense=1.5;
marks=[];
vcm=[];
type='';
thick=-Ngrid   ;
normfactor=1;  % option to enter a normalization for the data
titletag='';
labeltype='full';
titletype='full';
contSpacing=10;
contrCube=[];
contrType=[];
contLevs=[];
contColor='k';
contThick=thick;
streamColor='k';
vfieldColor='w';
markColor='k';
stlw=1;
slTypeDef='avg'; % is slice an average of slices or sum of them
contrWeight=ones(Ngrid,Ngrid,Ngrid);
xticks='';
yticks='';

axesHandle=0;
boxSize=0;
figureColorBk=[0 0 0];
figureColorText=[1 1 1];
newStarColor=[0.09,0.63,0.86];
lengthUnit='kpc';
guassSigma=0.5;

%Flags
cubeFlag=false;
typeFlag=false;
contrFlag=false;
wtFlag=false;
%labelFlag=false;
titleFlag=true;
printFlag=false;
clim_Flag=false; % no colorbar limits given, use min max
velFlag=false;
strmFlag=false;
%mvirboxFlag=false;
logFlag=false;
comoveFlag=true;
vcubeFlag(1:3)=false;
rogFlag=false;
zoomFlag=false;
contShowLabel=false;
contrLogFlag=false;
gridFlag=true;
rewerFlag=false;
starStructFlag=false;
newFigFlag=true;
smoothFlag=false;
xLabFlag=true;
yLabFlag=true;
saveFigFlag=false;
newStarFlag=false;

pdfFlag='pdf';
imageSmoothFlag='none';
brewerFlag=false;

map = brewermap(256,'*Greys');


%% reading arguments

i=1;
while i<=length(varargin)
    switch lower(varargin{i})
        
        % arguments for calculating cube based on gas Structure
        case{'star','starsstruct','stars'}
            starStructFlag=true;
            i=i+1;
            starStruct=varargin{i};
        case{'datatype','type'} %prepare data cube of this type
            i=i+1;
            typeFlag=true;
            type=varargin{i};
        case{'ngrid','ng'}  % set the grid size for the uniform cube
            i=i+1;
            Ngrid=varargin{i};
        case{'lengthunt','hublength','/h'}
            lengthUnit=[lengthUnit '/h'];
        case{'mask','starmask'}
            i=i+1;
            starMask=varargin{i};
            % arguments for plotting a given cube
        case{'box','boxsize'}
            i=i+1;
            boxSize=varargin{i};
            % arguments for plotting a given cube
            
        case{'data','datacube','cubestructure','cube'} %plot a given cube;
            i=i+1;
            cubeFlag=true;
            cubeStr=varargin{i};
            cube=cubeStr.cube;
            printTypeTag='data';
        case{'weightcube','wtcube','wtdatacube','wt'} % weight for cube
            i=i+1;
            wtFlag=true;
            weight=varargin{i};
            
            % set whether to average or sum up slice
        case{'avg','average'}
            slType='avg';
        case{'sum'}
            slType='sum';
            
        case{'smooth','nnsmooth','nn'}
            i=i+1;
            smoothFlag=true;
            nn=varargin{i};
        case{'imagesmooth','imagefilter','filter'}
            imageSmoothFlag='filter';
            i=i+1;
            imageFilter=varargin{i};
            
        case{'gausssmooth','gauss'}
            imageSmoothFlag='gauss';
            i=i+1;
            gaussSigma=varargin{i};
            % arguments for overlaying contours
        case {'newstars','shownew'}
            newStarFlag=true;
            i=i+1;
            newStarThresh=varargin{i};
        case {'newstarscolor','newcolor'}
                i=i+1;
            newStarColor=varargin{i};
        case {'contour','contourtype'}
            i=i+1;
            contrType=varargin{i};
            contrFlag=true;
        case {'contourdata','contourcube'}
            i=i+1;
            contrFlag=true;
            contrCube=varargin{i};
            %contrcubeFlag=true;
        case 'contlog'
            contrLogFlag=true;
        case {'contvals','contlevels'}
            i=i+1;
            contLevs=varargin{i};
        case 'contlabels'
            i=i+1;
            answer=varargin{i};
            contShowLabel=strcmp(answer,'on') ;
        case 'contcolor'
            i=i+1;
            contColor=varargin{i};
        case 'contthick'
            i=i+1;
            contThick=varargin{i};
        case 'avgcont'
            slTypeCont='avg';
        case 'sumcont'
            slTypeCont='sum';
            
            % velocity field stuff
        case 'streamcolor'
            i=i+1;
            streamColor=varargin{i};
        case {'vfieldcolor','velcolor'}
            i=i+1;
            vfieldColor=varargin{i};
        case{'rocube','rog'} % rps cube essential for calculating contours
            i=i+1;
            rogFlag=true;
            rog=varargin{i};
        case{'vxcube'} % option to enter a pre-loaded velocity cube
            i=i+1;
            vcubeFlag(1)=true;
            vx=varargin{i};
        case{'vycube'}
            i=i+1;
            vcubeFlag(2)=true;
            vy=varargin{i};
        case{'vzcube'}
            i=i+1;
            vcubeFlag(3)=true;
            vz=varargin{i};
        case {'velflag','velocity','vfield'}
            velFlag=true;
            strmFlag=true;
        case {'dilute','arrowspace'}
            i=i+1;
            dilute=varargin{i};
        case {'streamdense','streamlines'}
            i=i+1;
            strmFlag=true;
            streamDense=varargin{i};
            if ischar(streamDense)
                strmFlag=~strcmp(streamDense,'off');
            elseif isnumeric(streamDense)
                strmFlag=streamDense>0;
            end
        case {'vcm'}
            i=i+1;
            vcm=varargin{i};
            
            % basic arguments for plotting
        case{'log','log10'}
            logFlag=true;
        case{'nolog','linear'}
            logFlag=false;
        case{'width','thick'} %thicknes of slice in ckpc/h comoving (like box)
            i=i+1;
            thick=varargin{i};
        case {'yz','zy'}
            plotproj(1)=1;
        case {'zx','xz'}
            plotproj(2)=1;
        case {'xy','yx'}
            plotproj(3)=1;
        case 'all'
            plotproj=[1 1 1];
        case {'zoom','zoombox'}
            i=i+1;
            zoomBox=varargin{i};
            zoomFlag=true;
            if ischar(zoomBox)
                zoomFlag=~strcmp(zoomBox,'off');
            end
        case {'clims','limits'}
            i=i+1;
            clim_Flag=true;
            clims=varargin{i};
            if diff(clims)==0
                clim_Flag=false;
            end
            
            %         case {'alfa','alpha'}
            %             i=i+1;
            %             alfa=varargin{i};
        case {'normalize','norm','factor','normfactor'}
            i=i+1;
            normfactor=varargin{i};
            
            % arguments for overlying other things
        case {'circles','drawcircle','circle','circ'}
            i=i+1;
            circStruct=varargin{i};
            if ~isstruct(circStruct)
                error('MKMAP - draw circle input must be a structure')
            end
            
        case {'grid'}
            gridFlag=true;
        case {'nogrid'}
            gridFlag=false;
            
        case{'blackfigure','black'}
            figureColorBk='k';
            figureColorText='w';
            
        case{'whitefigure','white'}
            figureColorBk='w';
            figureColorText='k';
        case 'drawline'
            i=i+1;
            linStruct=varargin{i};
            if ~isstruct(linStruct)
                error('MKMAP - draw line input must be a structure')
            end
        case 'arrow'
            i=i+1;
            arrowStruct=varargin{i};
            if ~isstruct(arrowStruct)
                error('MKMAP - arrow input must be a structure')
            end
            
            
        case 'text'
            i=i+1;
            textBoxStruct=varargin{i};
            if ~isstruct(textBoxStruct)
                error('MKMAP - text box input must be a structure')
            end
        case {'brewer','brewermap'}
            i=i+1;
            brewMap=varargin{i};
            brewerFlag=true;
        case {'cmap','colormap'}
            i=i+1;
            map=varargin{i};
            brewerFlag=false;
            if size(map,2)~=3
                error('mkmapStars: Illegal colormap')
            end
            
        case 'xticks'
            i=i+1;
            xticks=varargin{i};
        case 'yticks'
            i=i+1;
            yticks=varargin{i};
        case 'fullticks'
            i=i+1;
            tickStruct=varargin{i};
            if ~isstruct(tickStruct)
                error('MKMAP - fullTick input must be a structure')
            end
        case {'noticks','notick'}
            xticks=' ';
            yticks=' ';
        case {'labels'}
            i=i+1;
            labeltype=varargin{i};
            %labelFlag=(strcmp(l,'yes') || strcmp(answer,'on') || strcmp(answer));
        case {'xlaboff','noxlab','noxlabel','xlabeloff'}
            xLabFlag=false;
        case {'ylaboff','noylab','noylabel','ylabeloff'}
            yLabFlag=false;
        case {'title'}
            i=i+1;
            titletag=varargin{i};
            titleFlag=~(strcmp(titletag,'no') || strcmp(titletag,'off') || strcmp(titletag,'none'));
        case {'titletype'}
            i=i+1;
            titletype=varargin{i};
        case {'bartag','bartitle'}
            i=i+1;
            bartagOveride=varargin{i};
             case {'barproperties' 'barprop'}
            i=i+1;
            barPropStruct=varargin{i};
            %          case {'mvirbox'} % place a text of Mvir on the map
            %             mvirboxFlag=true;
            %             bartag=varargin{i};
            %          case {'mvirbox'} % place a text of Mvir on the map
            %             mvirboxFlag=true;
            %             strpos=varargin{i};  % position of tex            strpos=varargin{i};  % position of text in map coordinates
        case {'comoving','comove','comoving_coord'}
            comoveFlag=true;
        case {'proper','proper_coord'}
            comoveFlag=false;
            
        case{'handle','fighandle','fig','figure'}
            i=i+1;
            hf=varargin{i};
            newFigFlag=false;
        case{'axes','axeshandle'}
            i=i+1;
            axesHandle=varargin{i};
            newFigFlag=false;
            
        case {'print'}
            printFlag=true;
            i=i+1;
            printtag=varargin{i};
        case {'nopdf','png'}
            pdfFlag='nopdf';
        case {'outputdir','printout','printoutdir'}
            i=i+1;
            printoutdir=varargin{i};
        case {'savefig'}
            %             i=i+1;
            %             printtag=varargin{i};
            saveFigFlag=true;
        otherwise
            error('mkmapStars: illegal argument %s',varargin{i})
    end
    i=i+1;
end


if ~any(plotproj)
    plotproj=[1 1 1];
end

if ~cubeFlag && ~(starStructFlag && typeFlag)
    error('mkmapStars: must enter both gas struct and data type to plot');
end

if ~(cubeFlag || typeFlag)
    error('mkmapStars: must enter data or datatype')
end
if (cubeFlag && typeFlag)
    error('mkmapStars: too many data arguments')
end

if contrFlag
    if ( isempty(contrCube) && isempty(contrType))
        error('mkmapStars: must enter contoure data or type')
    end
    
    if (~isempty(contrCube) && ~isempty(contrType))
        error('mkmapStars: too many contour data arguments')
    end
end

if ~wtFlag
    weight=ones(Ngrid,Ngrid,Ngrid);
end




%% enforce mask
noWindMask=starStruct.GFM_StellarFormationTime>0;
if isempty(starMask)
    starMask=noWindMask; %true(1,starStruct.count);
else
    starMask=starMask & noWindMask;
end

coord=starStruct.newCoord(:,starMask);
mass=starStruct.Masses(starMask);
%cellCount=sum(starMask);

if smoothFlag
    [~,dd]=knnsearch(coord',coord','K',nn+1);
    
    cellSize=dd(:,end)';
    
else
    cellSize=0;%(mass./starStruct.Density(starMask)).^(1/3);
end


%% build data cube
if typeFlag
    switch lower(type)
        % case {'flux','lightflux','light'}
        
        case {'density','surfacedensity','mass'}
            
            if sum(starMask)~=0
                
                cubeStr=cell2grid(coord,mass,cellSize,'ngrid',Ngrid,'extensive','box',boxSize);
                boxSize=cubeStr.boxSide;
                
                %cube=(cubeStr.cube.*massUnit)./(cubeStr.cellVol.*lengthUnit^3);
                cube=(cubeStr.cube)./(cubeStr.boxSide./cubeStr.Ngrid).^2.*illUnits.surfaceDensityUnit; % in Msun/kpc^2
                
            else
                cube=zeros(Ngrid,Ngrid,Ngrid);
            end
            
            weight=ones(size(cube));
            logFlag=true;
            bartag='$\log \rho_\mathrm{star}\,[\mathrm{M_\odot/kpc^2}]$';
            slTypeDef='sum';
            
            printTypeTag='dens';
            
            
        otherwise
            error('mkmapStars - unknown data type: %s',type)
    end
end
cube=cube.*normfactor;

%% find new star positions 
if newStarFlag
        nsZred=double(1./starStruct.GFM_StellarFormationTime(starMask)-1); 
        nsBirthTime=redshift2time(nsZred,'cosmo',cosmoStruct);  
        nsBirthTime=nsBirthTime.age;    
        currentAge=redshift2time(illUnits.zred,'cosmo',cosmoStruct);
        currentAge=currentAge.age;
        deltaAge=currentAge-nsBirthTime;
        if any(deltaAge<-1e-3); error('mkmapStars: something funky with the ages - probably mask is incorrect re wind particles');end
            
        newStarMask=deltaAge<=newStarThresh;
        newStarPos=coord(:,newStarMask);
end

%% set up contour data cube
% if ~isempty(contrType)
%     switch lower(contrType)
%         case{'none','off'}
%             contrFlag=false;
%         case {'density','gasdensity','rho','gas','rhogas'}
%             contrCube=rog;
%             contrWeight=ones(size(contrCube));
%             contrLogFlag=true;
%         case {'dmdensity','rhodm','dm'}
%             contrCube=smooth3(RHODM(boxx),'gaussian',9,5);
%
%             contrWeight=ones(size(contrCube));
%             contrLogFlag=true;
%         case {'totaldensity','rhotot','tot'}
%             contrCube=RHOTOT(boxx);
%             contrWeight=arones(size(contrCube));
%             contrLogFlag=true;
%         case {'entropy','ent','k','s'}
%             %SVIR=TVIR.*RVIR.^2./(3.*MVIR./(4.*pi)).^(2./3); %normalization for entropy
%             contrCube=S(boxx).*f_ent; % in units of KeV cm^2
%             contrWeight=rog;
%             contrLogFlag=true;
%         case {'temperature','temp','t'}
%             contrCube=T(boxx);
%             contrWeight=rog;
%             contrLogFlag=true;
%
%         case {'pressure','press','p'}
%             contrCube=T(boxx).*rog./(TVIR*MVIR/(4*pi*RVIR^3/3));
%             contrLogFlag=true;
%             contrWeight=ones(size(contrCube));
%
%         case {'xray'}
%
%             %tb=T(boxx);
%             %tb(tb<1e5)=0;
%             vol=(boxx/NCELL/hub/(1+zred))^3;
%             contrCube=cooling_cube(boxx,'cgs/Mpc').*vol; %rog.^2.*sqrt(tb);
%             contrWeight=ones(size(cube));
%             contrLogFlag=true; %cube=log10(cube);
%             slTypeCont='sum';
%         case {'metals','zia'}
%             contrCube=ZIa(boxx);
%
%         case {'iametals','zii'}
%             contrCube=ZII(boxx);
%
%         case {'metalicity','iimetals','ztot'}
%             contrCube=ZIa(boxx)+ZII(boxx);
%
%         otherwise
%             error('MKMAP - Illegal contour tpye: %s',contrType);
%     end
% end


%define slice index
if thick>0
    thk=ceil(0.5.*thick./boxSize.*Ngrid);   %% thick is in comoving Mpc/h
else
    thk=0.5*abs(thick);   %% defualt value of 6 cells for slice
end

% if contThick>0
%     cthk=ceil(0.5.*contThick./cubeStr.boxSide.*cubeStr.Nrgid);   %% thick is in comoving Mpc/h
% else
%     cthk=thk;
% end

slind=(floor(0.5.*Ngrid)-thk+1):1:(ceil(0.5.*Ngrid)+thk);
%slindCont=(floor(0.5.*Ngrid)-cthk+1):1:(ceil(0.5.*Ngrid)+cthk);

sideInd=1:Ngrid;
halfSide=0.5*boxSize;

side=[-halfSide halfSide];

% if(~comoveFlag)
%     side=side./(1+zred);
%     halfSide=halfSide./(1+zred);
% end

if ~exist('slType','var') % choose whether to average or sum slice
    slType=slTypeDef;
end

%% perpare slice to be plotted.
slice=mk_slice(cube,weight,slind,slType);
if logFlag
    slice=log10(slice);
end


%% prepare contours
% if ~exist('slTypeCont','var')
%     slTypeCont=slTypeDef;
% end
% % prepare contour slice
% if contrFlag
%     contrSlice=mk_slice(contrCube,contrWeight,slindCont,slTypeCont);
%     if contrLogFlag
%         contrSlice=log10(contrSlice);
%     end
% end
% cBsize=diff(side)./NCELL;
% cside=side(1)+0.5*cBsize:cBsize:side(2)-0.5*cBsize;


%% velocity field parameters (default filte=4)
% if velFlag2
%     [v_x, v_y]=mk_vfield(vx,vy,vz,rog,slind);
%     clear vx vy vz;
%     diluted_len = length(1:dilute:Ngrid);
%     diluted_jump = 2.*halfSide/(diluted_len-1);
%     notdiluted_jump = 2.*halfSide/(Ngrid-1);
%     [xxv, yyv] = meshgrid(-halfSide:diluted_jump:halfSide, -halfSide:diluted_jump:halfSide);
%     [xxs, yys] = meshgrid(-halfSide:notdiluted_jump:halfSide,-halfSide:notdiluted_jump:halfSide);
% end



%% draw map
if ~clim_Flag
    clims(1) = min(slice(slice>-Inf));
    clims(2) = max(slice(slice<Inf));
    
end
for projection = 1:3
    if plotproj(projection)
        
        
        if newFigFlag
            if ~exist('hf')
                hf=figure;
            else
                hf(end+1)=figure;
            end
            set(gcf,'color',figureColorBk);
        elseif exist('hf')
            figure(hf)
        else
            error('MKMAPGAS: No valid Figure handle given');
        end
        
        
        
        if axesHandle~=0
            axes(axesHandle);
        end
        
        % fix to avoid problems with setting the colorbar
        %         if logFlag
        %             clims(isinf(clims))=-30;
        %         end
        
        %% start plotting
        
        imj=squeeze(slice(:,:,projection));
        
        switch imageSmoothFlag
            case 'filter'
                imj=imfilter(imj,imageFilter);
            case 'gauss'
                imj=imgaussfilt(imj,gaussSigma);
        end
        
        
        imagesc(side,side,imj,clims);%axis equal tight;
        hold on
        
       
        
       
        %         %% add velocity field
        %         if velFlag2
        %             resampled_v_x = imresize(squeeze(v_x(:,:,projection)), [diluted_len diluted_len], 'bilinear');
        %             resampled_v_y = imresize(squeeze(v_y(:,:,projection)), [diluted_len diluted_len], 'bilinear');
        %             if strmFlag
        %                 h = streamslice(xxs,yys,squeeze(v_x(sideInd,sideInd,projection)),squeeze(v_y(sideInd,sideInd,projection)),streamDense);
        %                 set(h,'color',streamColor)
        %             end
        %             if velFlag
        %                 quiver(xxv,yyv,resampled_v_x,resampled_v_y, 1.3,'color',vfieldColor);
        %             end
        %         end
        
        %% plot new stars
        if newStarFlag
            switch projectiontags{projection}
                case 'YZ'
                    %plot(newStarPos(2,:),newStarPos(3,:),'.','color',newStarColor);
                    scatter(newStarPos(3,:),newStarPos(2,:),5,'markerfacecolor',newStarColor,'markerfacealpha',0.25,'markeredgealpha',0.0);
                case 'XZ'
                    %plot(newStarPos(1,:),newStarPos(3,:),'.','color',newStarColor);
                    scatter(newStarPos(1,:),newStarPos(3,:),5,'markerfacecolor',newStarColor,'markerfacealpha',0.25,'markeredgealpha',0.0);
                case 'XY'
                    scatter(newStarPos(1,:),newStarPos(2,:),5,'markerfacecolor',newStarColor,'markerfacealpha',0.25,'markeredgealpha',0.0);
                    %plot(newStarPos(1,:),newStarPos(2,:),'.','color',newStarColor);
            end
        end
        
        %% add contours
        
        %         if contrFlag
        %             scontrSlice=squeeze(contrSlice(:,:,projection));
        %
        %             if isempty(contLevs)
        %                 contLevs=linspace(min(min(scontrSlice)),max(max(scontrSlice)),contSpacing);
        %             elseif length(contLevs)==1
        %                 contLevs=linspace(min(min(scontrSlice)),max(max(scontrSlice)),contLevs);
        %             end
        %
        %             [~,ch]=contour(cside,cside,scontrSlice,contLevs,contColor,'linewidth',2);
        %             if contShowLabel
        %                 set(ch,'ShowText','on','TextStep',get(ch,'LevelStep')*2)
        %             end
        %             %clabel(C,chh,'fontsize',18,'rotation',0)
        %         end
        
        %% draw circles
        if exist('circStruct','var')
            for indC=1:length(circStruct) % go over all lines
                
                if ~isfield(circStruct(indC),'color')||isempty(circStruct(indC).color) % set color  if none is given
                    circStruct(indC).color=figureColorText;
                end
                if ~isfield(circStruct(indC),'width')||isempty(circStruct(indC).width) % set linewidth  if none is given
                    circStruct(indC).width=1.5;
                end
                if ~isfield(circStruct(indC),'type')||isempty(circStruct(indC).type) % set line type  if none is given
                    circStruct(indC).type='-';
                end
                if ~isfield(circStruct(indC),'center')||isempty(circStruct(indC).center) % set line type  if none is given
                    circStruct(indC).center=[0 0];
                end
                
                
                drawCircle(hf(end),circStruct(indC).radius,...
                    'center',circStruct(indC).center,...
                    'clr',circStruct(indC).color,...
                    'lw',circStruct(indC).width,...
                    'style',circStruct(indC).type)
            end
            
        end
        hold off
        
        
        
        
        
        
        %         % make sure marked circles fit in box
        %         if ~isempty(marks)
        %             if(comoveFlag)
        %                 markss=marks.*(1+zred);
        %             else
        %                 markss=marks;
        %             end
        %
        %             drawCircles(hf,marks_inbox,markColor);
        %         end
        %
        %         %                 rv(1)=get_rvir(500);
        %         %                 rv(2)=get_rvir();
        %         %                 rv(3)=get_rvir(200);
        %         %
        %         %                 if(comoveFlag)
        %         %                     rv=rv.*(1+zred);
        %         %                 end
        %         %                 draw_circle_boxx(gcf,rv(rvin),'white');      % draw rvir
        
        
        
        
        
        
        %% draw a line
        if exist('linStruct','var')
            for indL=1:length(linStruct) % go over all lines
                
                if ~isfield(linStruct(indL),'color')||isempty(linStruct(indL).color) % set color  if none is given
                    linStruct(indL).color='k';
                end
                if ~isfield(linStruct(indL),'width')||isempty(linStruct(indL).width) % set linewidth  if none is given
                    linStruct(indL).width=1.5;
                end
                if ~isfield(linStruct(indL),'type')||isempty(linStruct(indL).type) % set line type  if none is given
                    linStruct(indL).type='-';
                end
                linVal=linStruct(indL).value;
                if length(linVal)==1
                    linVal=linVal.*[1 1];
                end
                %draw line
                switch(linStruct(indL).dir)
                    case 'horizontal'
                        xl=xlim;
                        yl=linVal;
                        
                    case 'vertical'
                        yl=ylim;
                        xl=linVal;
                        
                    otherwise
                        error('MKMAP - Illegal direction in draw line: %s',linStruct.dir)
                end
                hold on
                plot(xl,yl,linStruct(indL).type,'color',linStruct(indL).color,'linewidth',linStruct(indL).width);
                hold off
            end
            
        end
        hold off
        
        %% draw an arrow
        if exist('arrowStruct','var')
            for indA=1:length(arrowStruct) % go over all lines
                if ~isfield(arrowStruct(indA),'start')  % set line type  if none is given
                    error('MKMAP: Arrow structure must have a START field')
                end
                if ~isfield(arrowStruct(indA),'stop') % set line type  if none is given
                    error('MKMAP: Arrow structure must have a STOP field')
                end
                if ~isfield(arrowStruct(indA),'color') || isempty(arrowStruct(indA).color) % set color  if none is given
                    arrowStruct(indA).color='k';
                end
                if ~isfield(arrowStruct(indA),'width') || isempty(arrowStruct(indA).width)% set linewidth  if none is given
                    arrowStruct(indA).width=1.5;
                end
                if ~isfield(arrowStruct(indA),'faceColor') || isempty(arrowStruct(indA).faceColor) % set linewidth  if none is given
                    arrowStruct(indA).faceColor=arrowStruct(indA).color;
                end
                if ~isfield(arrowStruct(indA),'edgeColor') || isempty(arrowStruct(indA).edgeColor)% set linewidth  if none is given
                    arrowStruct(indA).edgeColor=arrowStruct(indA).color;
                end
                if ~isfield(arrowStruct(indA),'ends') || isempty(arrowStruct(indA).ends)% set linewidth  if none is given
                    arrowStruct(indA).ends='stop';
                end
                hold on
                %                 start=arrowStruct(indA).start;
                %                 stop=arrowStruct(indA).stop;
                %                 xl=xlim;
                %                 yl=ylim;
                %                 xx=([start(1) stop(1)]-xl(1))./diff(xl);
                %                 yy=([start(2) stop(2)]-yl(1))./diff(yl);
                %
                %                 ha=annotation('arrow',xx,yy,...
                %                     'Color',arrowStruct(indA).faceColor,'Linewidth',arrowStruct(indA).width);
                %
                %
                arrow(arrowStruct(indA).start,arrowStruct(indA).stop,...
                    'width',arrowStruct(indA).width',...
                    'FaceColor',arrowStruct(indA).faceColor,...
                    'EdgeColor',arrowStruct(indA).edgeColor,...
                    'Ends',arrowStruct(indA).ends);
            end
        end
        
        hold off
        
        
        
        %% put axis labels
        prjtag=projectiontags{projection};
        
        if(comoveFlag && illUnits.zred>0)
            lunit=['[\mathrm{' lengthUnit '\, comoving\,}]'];
        else
            lunit=['[\mathrm{' lengthUnit '}]'];
        end
        
        
        
        switch labeltype
            case {'yes','on','full'}
                switch projection
                    case 1  % YZ
                        if xLabFlag
                            xlh=xlabelmine(sprintf('$Z %s $',lunit));
                            set(xlh,'color',figureColorText);
                        end
                        if yLabFlag
                            ylh=ylabelmine(sprintf('$Y %s $',lunit));
                            set(ylh,'color',figureColorText);
                        end
                    case 2  % ZX
                        if xLabFlag
                            xlh=xlabelmine(sprintf('$X %s $',lunit));
                            set(xlh,'color',figureColorText);
                        end
                        if yLabFlag
                            ylh=ylabelmine(sprintf('$Z %s $',lunit));
                            set(ylh,'color',figureColorText)
                        end
                    case 3  % XY
                        if xLabFlag
                            xlh=xlabelmine(sprintf('$X %s $',lunit));
                            set(xlh,'color',figureColorText);
                        end
                        if yLabFlag
                            ylh=ylabelmine(sprintf('$Y %s $',lunit));
                            set(ylh,'color',figureColorText);
                        end
                end
                
            case {'half','units'}
                if xLabFlag
                    xlh=xlabelmine(sprintf('$ %s $',lunit));
                    set(xlh,'color',figureColorText);
                end
                
                if yLabFlag
                    ylh=ylabelmine(sprintf('$ %s $',lunit));
                    set(ylh,'color',figureColorText);
                end
                
        end
        
        
        
        
        
        %% sort out ticks
        %         if ~exist('tickStruct','var')
        %             tickjump=0.25.*boxx;
        %
        %             if ~comoveFlag && strcmp(aexpn,'a06')
        %                 tickjump=0.15.*boxx;
        %             end
        %             ticks=-2*tickjump:tickjump:2*tickjump;
        %             if ~exist('xticks','var')
        %                 xticks=ticks;
        %             end
        %             if ~exist('yticks','var')
        %                 yticks=ticks;
        %             end
        %             if boxx==1
        %                 xtickLabels=sprintf('%3.2f|',xticks);
        %                 ytickLabels=sprintf('%3.2f|',yticks);
        %                 zl='|0.00|';
        %
        %             else  %if any(boxx==[2 4])
        %                 xtickLabels=sprintf('%3.1f|',xticks);
        %                 ytickLabels=sprintf('%3.1f|',yticks);
        %                 zl='|0.0|';
        %                 %             else
        %                 %                 xtickLabels=sprintf('%3.0f|',xticks);
        %                 %                 ytickLabels=sprintf('%3.0f|',yticks);
        %                 %                 zl='|0|';
        %                 %
        %             end
        %
        %             % fix zero in labels
        %             k=strfind(xtickLabels,zl);
        %             if k>0
        %                 xtickLabels=sprintf('%s%s%s',xtickLabels(1:k-1),'|0|',xtickLabels(k+length(zl):end));
        %             end
        %
        %             k=strfind(ytickLabels,zl);
        %             if k>0
        %                 ytickLabels=sprintf('%s%s%s',ytickLabels(1:k-1),'|0|',ytickLabels(k+length(zl):end));
        %             end
        %
        %         else
        %             xticks=tickStruct.xticks;
        %             yticks=tickStruct.yticks;
        %             %xtickLabels=tickStruct.xtickLabels;
        %             %ytickLabels=tickStruct.ytickLabels;
        %         end
        
        
        % 'XTick',xticks ,'YTick',yticks,... %'TickLength',[-0.015 -0.015],...
        set(gca,'Ydir','normal','Fontsize',14,...
            'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1],...
            'XLim',side,'YLim',side,'color',figureColorText,...
            'xcolor',figureColorText,'ycolor',figureColorText);    %'XTickLabel',xtickLabels,'YTickLabel',ytickLabels
        
        if ~isempty(xticks)
            set(gca,'XTickLabel',xticks);
        end
        if ~isempty(yticks)
            set(gca,'YTickLabel',yticks);
        end
        
        box on;
        
        %% address Zoom
        if zoomFlag
            zb=[-1 1;-1 1];
            if length(zoomBox)==1
                zb=zb.*zoomBox;
            elseif length(zoomBox)==3
                zb(1,:)= zoomBox(1).*[1 1]+zoomBox(3).*[0 1];
                zb(2,:)= zoomBox(2).*[1 1]+zoomBox(3).*[0 1];
            else
                error('MKMAP - Illegal zoom box')
            end
            set(gca,'XLim',zb(1,:),'YLim',zb(2,:))
            
            
        end
        
        %% do bar and colormap stuff
        caxis(clims);
        
        if brewerFlag
            map = brewermap(256,brewMap);
        end
        %map(1,:)=[0 0 0];
        
        
        colormap(map);
        bar=colorbar;
        if ~isempty(bartagOveride) || ~exist('bartag','var')
            bartag=bartagOveride;
        end
        barTitle(bar,bartag,'color',figureColorText);

        if exist('barPropStruct','var')
            setPropertiesStructure(bar,barPropStruct)
        else
            set(bar,'Fontsize',14,'color',figureColorText)
        end
        
        
        
        
        %set(get(bar,'Title'),'String',bartag,'Fontsize',12,'Interpreter','latex');
        
        set(gcf,'Colormap',map);
        
        %% set title
        if titleFlag
            ht=titlemine(titletag);
            set(ht,'color',figureColorText)
        end
        %title(sprintf('%s %s, Thickness=%s Mpc/h',CLUSTER,type,num2str(thick,3)),'Fontsize',12,'Interpreter','latex');
        
        %         %% set title
        %         if titleFlag
        %             switch titletype
        %                 case 'full'
        %                     mvstr=mk_mvir_string(MVIR);
        %                     titlemine(sprintf('%s %s, %s $z=%3.2g$',CLUSTER,mvstr,titletag,zred));
        %                 case 'custom'
        %                     titlemine(titletag);
        %                 otherwise
        %                     error('MKMAP: Illegal value for title type (full/custom)')
        %             end
        %         end
        
        %% add grid
        if gridFlag
            grid;
        end
        
        
        %% add text box
        
        if exist('textBoxStruct','var')
%             axPos=get(gca,'pos');
%             xl=xlim;
%             yl=ylim;
            
            
            
            for kk=1:length(textBoxStruct)
                if ~isfield(textBoxStruct(kk),'fontsize') || isempty(textBoxStruct(kk).fontsize)% set linewidth  if none is given
                    textBoxStruct(kk).fontsize=20;
                end
                
                % calculate position of text (given in data units) to position
                % within the axes and figure, as needed by annotation (0-1 in
                % figure)
                pos=double(textBoxStruct(kk).position);

                text(pos(1),pos(2),textBoxStruct(kk).str,...
                    'color',textBoxStruct(kk).color,...
                    'Fontsize',textBoxStruct(kk).fontsize,...
                     'Interpreter','latex')
                
                
%                 tpos(1)=axPos(1)+axPos(3)*(pos(1)-xl(1))./diff(xl);
%                 tpos(2)=axPos(2)+axPos(4)*(pos(2)-yl(1))./diff(yl);
%                 tpos(3)=axPos(3)*pos(3)./diff(xl);
%                 tpos(4)=axPos(4)*pos(4)./diff(yl);
%                 
%                 %textBoxStruct.pos=
%                 annotation(gcf,'textbox',double(tpos),...
%                     'Color',textBoxStruct(kk).color,...
%                     'String', textBoxStruct(kk).str,...
%                     'Interpreter','latex',...
%                     'Fontsize',textBoxStruct(kk).fontsize,...
%                     'FitBoxToText','off','linestyle','none');
                
            end
        end
        
        %% print
        if printFlag
            
            name=sprintf('starMap%s_%s_%s_snp%i_%s',prjtag,printTypeTag,printtag,illUnits.snap,simDisplayName);
            %name=sprintf('%s/%s_map%s_b%d_%s.%s',printoutdir,CLUSTER,prjtag,boxx,printtag,'%s');
            
            
            printout_fig(gcf,name,'dir',printoutdir,'v',pdfFlag)
            
            
        end
        
        
        
    end
end

if saveFigFlag
    
    name=sprintf('%s/figFiles/starMap_%s_%s_snp%i_%s.fig',DEFAULT_PRINTOUT_DIR,printTypeTag,printtag,illUnits.snap,simDisplayName);
    %name=sprintf('%s/%s_map%s_b%d_%s.%s',printoutdir,CLUSTER,prjtag,boxx,printtag,'%s');
    
    savefig(hf,name,'compact')
    
end


%
%

end

