function  hf=mkmapGas(varargin )
%MKMAP plotting the TNG objects
%   Plotting the gas cells from the TNG by mapping to a uniform grid
%


warning('%s - This function does not fix the length units - make sure your coordinates are scales accordingly!',...
    current_function().upper);

%% defuals and globals
global DEFAULT_PRINTOUT_DIR
global DEFAULT_MATFILE_DIR
global illUnits
global simDisplayName

units;

zred=0;
projectiontags={'YZ' 'XZ' 'XY'};
bartagOveride='';
printoutdir=DEFAULT_PRINTOUT_DIR;
printtag='badtag';
Ngrid=256;
plotproj=[0 0 0]; % plot all three projections
dilute=4;
streamDense=1.5;
marks=[];
vcm=[0,0,0];
type='';
thick=-Ngrid   ;
normfactor=1;  % option to enter a normalization for the data
titletag='';
labeltype='full';
%titletype='full';
contSpacing=10;
contrCube=[];
contrType=[];
contLevs=[];
contColor=[0 0 0];
contThick=thick;
streamColor=[0 0 0];
vfieldColor=[1 1 1];
markColor=[0 0 0];
stlw=1;
slTypeDef='avg'; % is slice an average of slices or sum of them
contrWeight=ones(Ngrid,Ngrid,Ngrid);
%hf=0;
axesHandle=0;
boxSize=0;
gasMask=[];
figureColorBk=[1 1 1]; %'w';
figureColorText=[ 0 0 0];%'k';
nanVal='none';
lengthUnit='kpc/h';
idObj=[];
r200c=[];
guassSigma=0.5;
xticks='';
yticks='';

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
brewerFlag=false;
gasStructFlag=false;
newFigFlag=true;
backgroundFlag=false; %true;
map = brewermap(256,'*Spectral');
saveFigFlag=false;
visibleFlag=true;
pointPlotFlag=false;
imageSmoothFlag='none';
xLabFlag=true;
yLabFlag=true;
pdfFlag='pdf';

%% reading arguments

i=1;
while i<=length(varargin)
    switch lower(varargin{i})
        
        % arguments for calculating cube based on gas Structure
        case{'gas','gasstruct'}
            gasStructFlag=true;
            i=i+1;
            gasStruct=varargin{i};
        case{'datatype','type'} %prepare data cube of this type
            i=i+1;
            typeFlag=true;
            type=varargin{i};
        case{'id','idgal','idfof','idobj','object'}
            i=i+1;
            idObj=varargin{i};
        case{'r200','r200c','r200crit'}
            i=i+1;
            r200c=varargin{i};
        case{'ngrid','ng'}  % set the grid size for the uniform cube
            i=i+1;
            Ngrid=varargin{i};
        case{'lengthunt','hublength','/h'}
            lengthUnit=[lengthUnit '/h'];
        case{'fixed length','no hub length','no /h','no\h'}
            lengthUnit='kpc';
        case{'lengthunit'}
            i=i+1;
            lengthUnit=varargin{i};
        case{'mask','gasmask'}
            i=i+1;
            gasMask=varargin{i};
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
        case{'imagesmooth','imagefilter','filter'}
            imageSmoothFlag='filter';
            i=i+1;
            imageFilter=varargin{i};
            
         case{'gausssmooth','gauss'}
            imageSmoothFlag='gauss';
            i=i+1;
            gaussSigma=varargin{i};   
            % arguments for overlaying contours
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
            if length(vcm)~=3
               error('%s - vcm must be a 3-component array',current_function().upper)
            end
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
        case {'nanval','nanvalues'}
            i=i+1;
            nanVal=varargin{i};
        case {'normalize','norm','factor','normfactor'}
            i=i+1;
            normfactor=varargin{i};
            
            % arguments for overlying other things
        case {'circles','drawcircle','circle','circ'}
            i=i+1;
            circStruct=varargin{i};
            if ~isstruct(circStruct)
                error('%s - draw circle input must be a structure',current_function().upper)
            end
            
        case {'grid'}
            gridFlag=true;
        case {'nogrid'}
            gridFlag=false;
        case{'blackfigure','black'}
            figureColorBk=[0 0 0];
            figureColorText=[1 1 1];
            
        case 'drawline'
            i=i+1;
            linStruct=varargin{i};
            if ~isstruct(linStruct)
                error('%s - draw line input must be a structure',current_function().upper)
            end
        case 'arrow'
            i=i+1;
            arrowStruct=varargin{i};
            if ~isstruct(arrowStruct)
                error('%s - arrow input must be a structure',current_function().upper)
            end
        case 'text'
            i=i+1;
            textBoxStruct=varargin{i};
            if ~isstruct(textBoxStruct)
                error('%s - text box input must be a structure',current_function().upper)
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
                error('%s - Illegal colormap',current_function().upper)
            end
        case {'nobackground'}
            backgroundFlag=false;
            case {'background'}
            backgroundFlag=true;
             case {'backgroundcolor','bkcolor','colorbk'}
            i=i+1;
            figureColorBk=varargin{i};
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
                error('%s - fullTick input must be a structure',current_function().upper)
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
            %         case {'titletype'}
            %             i=i+1;
            %             titletype=varargin{i};
        case {'bartag','bartitle'}
            i=i+1;
            bartagOveride=varargin{i};
        case {'barproperties' 'barprop'}
            i=i+1;
            barPropStruct=varargin{i};
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
            
            %case{'invisible','noshow','noplot'}
            
            %% printing
        case {'print'}
            printFlag=true;
            i=i+1;
            printtag=varargin{i};
        case {'outputdir','printout','printoutdir'}
            i=i+1;
            printoutdir=varargin{i};
        case {'savefig'}
            %             i=i+1;
            %             printtag=varargin{i};
            saveFigFlag=true;
        case {'nopdf','png'}
            pdfFlag='nopdf';
        otherwise
            error('%s - illegal argument: %s',current_function().upper,varargin{i})
    end
    i=i+1;
end


if ~any(plotproj)
    plotproj=[1 1 1];
end

if ~cubeFlag && ~(gasStructFlag && typeFlag)
    error('%s - must enter both gas struct and data type to plot',current_function().upper);
end

if ~(cubeFlag || typeFlag)
    error('%s - must enter data or datatype',current_function().upper)
end
if (cubeFlag && typeFlag)
    error('%s - too many data arguments',current_function().upper)
end

if contrFlag
    if ( isempty(contrCube) && isempty(contrType))
        error('%s - must enter contoure data or type',current_function().upper)
    end
    
    if (~isempty(contrCube) && ~isempty(contrType))
        error('%s - too many contour data arguments',current_function().upper)
    end
end

if ~wtFlag
    weight=ones(Ngrid,Ngrid,Ngrid);
end

%% add missing fields
if (gasStructFlag)
    
    if ~isfield(gasStruct,'newCoord')
       error('%s - object has not been centered !!',current_function().upper)
    end
end


%% enforce mask
if isempty(gasMask)
    gasMask=true(1,gasStruct.count);
end
if isfield(gasStruct,'StarFormationRate')
sfrMask=gasStruct.StarFormationRate(gasMask)==0;
end
coord=gasStruct.newCoord(:,gasMask);
mass=gasStruct.Masses(gasMask);
cellSize=2.*(3.*mass./(4.*pi.*gasStruct.Density(gasMask))).^(1/3); % "diameter" of a cell 
cellCount=sum(gasMask);




%% build data cube
if typeFlag
    switch lower(type)
        case {'flux','massflux','f'}
            %fluxnorm=(0.1117.*(MVIR./1e15).^0.15*(1+zred)^2.25)./(4.*pi); % normalization for flux
            
            % find vr
            vv=zeros(3,cellCount);
            dist=sqrt(sum(coord.^2,1));
            vr=zeros(3,cellCount);
            for k=1:3
                vv(k,:)=gasStruct.Velocities(k,gasMask)-vcm(k);
                vr(k,:)=vv(k,:).*coord(k,:)./dist;
                
            end
            vrr=sum(vr,1);
            
            cubeStr=cell2grid(coord,vrr.*mass,cellSize,'ngrid',Ngrid,'extensive','box',boxSize);
            cube=cubeStr.cube./(cubeStr.cellVol).*illUnits.densityUnit.*illustris.utils.velocityFactor(illUnits.snap,'gas');  % Msun/kpc^3*km/sec
            logFlag=false;
            %weight=cubeStr.weights; % cube of mass in each uniform grid cell
            bartag='$\rho v_r\,[\mathrm{M_\odot/kpc^3\,km/sec}]$';
            slTypeDef='sum';
            map = brewermap(256,'*RdBu');
            
            %             %find radial velocity component
            %             [meshY, meshX, meshZ] = meshgrid(1:size(vy,1), 1:size(vx,2), 1:size(vz,3));
            %             %convert to center origin coordinates
            %             meshX = meshX - (size(vx,1)+1)/2 -cm(1);
            %             meshY = meshY - (size(vy,2)+1)/2 -cm(2);
            %             meshZ = meshZ - (size(vz,3)+1)/2 -cm(3);
            %             % Fix Units (to be in Mpc)
            %             meshX = meshX * ((boxx/hub)/NCELL);
            %             meshY = meshY * ((boxx/hub)/NCELL);
            %             meshZ = meshZ * ((boxx/hub)/NCELL);
            %
            %             rcube2=meshX.^2+meshY.^2+meshZ.^2 ; % r^2 cube in Mpc
            %             vr=(vx.*meshX+vy.*meshY+vz.*meshZ)./sqrt(rcube2) ; %radial velocity
            %
            %             switch SIMTYPE
            %                 case 'CSF'
            %                     [Mg200 , ~, ~]=read_Mass_Profiles(RVIR);
            %                 case 'Adiabatic'
            %                     bx=2^(ceil(log(2*get_rvir*hub)/log(2)));
            %                     mg=RHOG(bx).*(bx./NCELL./hub).^3;
            %                     mgP=MAKE_CUM_PROFILE_FROM_CUBE(mg);
            %                     ind=ceil(get_rvir/(0.5*bx/hub/NCELL));
            %                     Mg200=mgP(ind);
            %             end
            %
            %             cube=rog.*vr.*rcube2./Mg200.*(km/Mpc*Gyr)./fluxnorm ;
            %             weight=ones(size(cube));
            %             bartag='$\frac{\dot{M}_{gas}}{d\Omega}\,[\frac{\dot{M}_{gas}}{d\Omega}|_\mathrm{vir} ]$';
            %             %bartag='$\frac{\dot{M}/M_{gas}}{d\Omega} [\frac{\dot M}{M_{gas}}\left|_{\mathrm{virial}}]$';
            %             %bartag='$\dot{M}_{gas}$';
            %
            %             printTypeTag='flux';
            %
            %             clear vr rcube2 meshX meshY meshZ
            
            
        case {'density','gasdensity','rho','gas','rhogas'}
            cubeStr=cell2grid(coord,mass,cellSize,'ngrid',Ngrid,'extensive','box',boxSize);
            
            %cube=(cubeStr.cube.*massUnit)./(cubeStr.cellVol.*lengthUnit^3);
            cube=(cubeStr.cube)./(cubeStr.cellVol).*illUnits.densityUnit; % in Msun/kpc^3
            weight=ones(size(cube));
            logFlag=true;
            bartag='$\log \rho_{gas}\,[\mathrm{M_\odot/kpc^3}]$';
            slTypeDef='avg';
            printTypeTag='dens';
			
		
            
        case {'numberdensity','n','ndensity'}
            cubeStr=cell2grid(coord,mass,cellSize,'ngrid',Ngrid,'extensive','box',boxSize);
            
            %cube=(cubeStr.cube.*massUnit)./(cubeStr.cellVol.*lengthUnit^3);
            
            %rhoFac=illUnits.densityUnit.*(Units.Ms/Units.kpc^3/Units.mm); %in cm^-3
            cube=(cubeStr.cube)./(cubeStr.cellVol).*illUnits.numberDensityFactor; % in cm^-3
            weight=ones(size(cube));
            logFlag=true;
            bartag='$\log n\,[\mathrm{cm^{-3}}]$';
            slTypeDef='avg';
            
            printTypeTag='nDens';
            
%         case {'coldensity','columndensity','rhocol','gascol','rhogascol'}
%             cubeStr=cell2grid(coord,mass,cellSize,'ngrid',Ngrid,'extensive','box',boxSize);
%             
%             %cube=(cubeStr.cube.*massUnit)./(cubeStr.cellVol.*lengthUnit^3);
%             cube=(cubeStr.cube)./(cubeStr.).*illUnits.densityUnit; % in Msun/kpc^3
%             weight=ones(size(cube));
%             logFlag=true;
%             bartag='$\log \rho_{gas}\,[\mathrm{M_\odot/kpc^3}]$';
%             slTypeDef='sum';
%             printTypeTag='dens';
			
            
            
        case {'hi','hidensity'}
            mhi=gas.mHi_BR(gasMask);
            
            cubeStr=cell2grid(coord,mhi,cellSize,'ngrid',Ngrid,'extensive','box',boxSize);
            
            %cube=(cubeStr.cube.*massUnit)./(cubeStr.cellVol.*lengthUnit^3);
            
            
            cube=(cubeStr.cube)./(cubeStr.cellVol).*illUnits.densityUnit.*... % this        is in Msun/kpc^3
                (illUnits.physUnits.Ms/illUnits.physUnits.kpc^3/illUnits.physUnits.mp); % and this changes to atoms/cm^3
            weight=ones(size(cube));
            logFlag=true;
            bartag='$\log n_\mathrm{HI}\,[\mathrm{cm^{-3}}]$';
            slTypeDef='avg';
            
            printTypeTag='hiDens';
        case {'hicol','hicoldensity'}
            mhi=gasStruct.mHi_BR(gasMask);
            
            cubeStr=cell2grid(coord,mhi,cellSize,'ngrid',Ngrid,'extensive','box',boxSize);
            
            cube=(cubeStr.cube)./(cubeStr.cellArea).*illUnits.surfaceDensityUnit.*... % this is in Msun/kpc^2
                (illUnits.physUnits.Ms/illUnits.physUnits.kpc^2/illUnits.physUnits.mp); % and this changes to atoms/cm^2
            weight=ones(size(cube));
            logFlag=true;
            bartag='$\log n_\mathrm{HI,col}\,[\mathrm{cm^{-2}}]$';
            slTypeDef='sum';
            
            printTypeTag='hiColDens';
            
        case {'h2','h2density'}
            mh2=gas.mH2_BR(gasMask);
            
            cubeStr=cell2grid(coord,mh2,cellSize,'ngrid',Ngrid,'extensive','box',boxSize);
            
            %cube=(cubeStr.cube.*massUnit)./(cubeStr.cellVol.*lengthUnit^3);
            
            
            cube=(cubeStr.cube)./(cubeStr.cellVol).*illUnits.densityUnit.*... % this        is in Msun/kpc^3
                (illUnits.physUnits.Ms/illUnits.physUnits.kpc^3/(2*illUnits.physUnits.mp)); % and this changes to atoms/cm^3
            weight=ones(size(cube));
            logFlag=true;
            bartag='$\log n_\mathrm{HI}\,[\mathrm{cm^{-3}}]$';
            slTypeDef='avg';
            
            printTypeTag='hiDens';
        case {'h2col','h2coldensity'}
            mh2=gasStruct.mH2_BR(gasMask);
            
            cubeStr=cell2grid(coord,mh2,cellSize,'ngrid',Ngrid,'extensive','box',boxSize);
            
            cube=(cubeStr.cube)./(cubeStr.cellArea).*illUnits.surfaceDensityUnit.*... % this is in Msun/kpc^2
                (illUnits.physUnits.Ms/illUnits.physUnits.kpc^2/(2*illUnits.physUnits.mp)); % and this changes to atoms/cm^2
            weight=ones(size(cube));
            logFlag=true;
            bartag='$\log n_\mathrm{HI,col}\,[\mathrm{cm^{-2}}]$';
            slTypeDef='sum';
            
            printTypeTag='hiColDens';
            
            
        case {'entropy','ent','k','s'}
            if ~isfield(gasStruct,'Entropy')
                gasStruct=illustris.utils.addEntropy(gasStruct); %in KeV cm^2
            end
            
            ent=gasStruct.Entropy(gasMask);
            msk=sfrMask;
            
            cubeStr=cell2grid(coord(:,msk),ent(msk),cellSize(msk),...
                'ngrid',Ngrid,'intensive','weights',mass(msk),'box',boxSize);
            cube=cubeStr.cube; %in K
            logFlag=true;
            weight=cubeStr.weights; % cube of mass in each uniform grid cell
            bartag='$ S\,[\mathrm{KeV\, cm^2}]$';
            slTypeDef='avg';
            printTypeTag='ent';
            
        case {'tcool','coolingtime','tc'}
            
            tcool=illustris.utils.calcCoolingTime(gasStruct.Density,gasStruct.InternalEnergy,gasStruct.GFM_CoolingRate); % cooling time in Gyr^-1
            
            % avoid cells with tc=0;
            tcool=tcool(gasMask);
            msk=(tcool>0 & sfrMask);
            
            cubeStr=cell2grid(coord(:,msk),tcool(msk),cellSize(msk),...
                'ngrid',Ngrid,'intensive','weights',mass(msk),'box',boxSize);
            cube=cubeStr.cube; %in K
            logFlag=true;
            weight=cubeStr.weights; % cube of mass in each uniform grid cell
            bartag= '$\log t_\mathrm{cool}\,[\mathrm{Gyr}]$' ;
            slTypeDef='avg';
            printTypeTag='tcool';
            
        case {'tctff','coolfreefall'}
            
            tcool=illustris.utils.calcCoolingTime(gasStruct.Density,gasStruct.InternalEnergy,gasStruct.GFM_CoolingRate); % cooling time in Gyr^-1
            
            % avoid cells with tc=0;
            tcool=tcool(gasMask);
            msk=(tcool>0 & sfrMask);
            
            cubeStr=cell2grid(coord(:,msk),tcool(msk),cellSize(msk),...
                'ngrid',Ngrid,'intensive','weights',mass(msk),'box',boxSize);
            cube=cubeStr.cube; %in K
            
            logFlag=true;
            weight=cubeStr.weights; % cube of mass in each uniform grid cell
            
            load([DEFAULT_MATFILE_DIR '/freeFallTime_profiles_snp' num2str(illUnits.snap) '_' simDisplayName '.mat'])
            if isempty(idObj)
               error('%s - no ID - needed for tff calcultion',current_function().upper);
            end
            
            if isempty(r200c)
               error('%s - no r200c - needed for tff calcultion',current_function().upper);
            end
            
            rc=generate_radius_cube('ng',cubeStr.Ngrid,'size',double(cubeStr.boxSide));
            aCof=tffProfile.polyfit.a(idObj+1);
            bCof=tffProfile.polyfit.b(idObj+1);
            cCof=tffProfile.polyfit.c(idObj+1);
            
            tffCube=recreate_TFF(rc,double(r200c),aCof,bCof,cCof);
            
            cube=cube./tffCube;
            
            bartag= '$\log t_\mathrm{cool}/t_\mathrm{ff}$' ;
            slTypeDef='avg';
            printTypeTag='tctff';
            
            clear tffCube rc
        case {'temperature','temp','t'}
            
            if ~isfield(gasStruct,'Temperature')
                gasStruct=illustris.utils.addTemperature(gasStruct); %in K
            end
            
            temp=gasStruct.Temperature(gasMask);
            msk=sfrMask;
            
            cubeStr=cell2grid(coord(:,msk), temp(msk),cellSize(msk),...
                'ngrid',Ngrid,'intensive','weights',mass(msk),'box',boxSize);
            cube=cubeStr.cube; %in K
            logFlag=true;
            weight=cubeStr.weights; % cube of mass in each uniform grid cell
            bartag='$\log T\,[\mathrm{K}]$';
            slTypeDef='avg';
            printTypeTag='temp';
            
        case {'mach'}
            
            cubeStr=cell2grid(coord,gasStruct.Machnumber(gasMask),ones(size(cellSize)),...
                'ngrid',Ngrid,'intensive','weights',gasStruct.EnergyDissipation(gasMask),'box',boxSize);
            cube=cubeStr.cube; %in K
            logFlag=true;
            weight=cubeStr.weights; % cube of mass in each uniform grid cell
            bartag='$\mathcal{M}$';
            slTypeDef='avg';
            printTypeTag='mach';
            
            map=brewermap(256,'*YlOrRd');
            map(1,:)=[0 0 0];
            backgroundFlag=false;
            
        case {'shock','shockenergy','energydissipation','dissipation'}
            cubeStr=cell2grid(coord,gasStruct.EnergyDissipation(gasMask),zeros(size(cellSize)),...
                'ngrid',Ngrid,'max','box',boxSize);
            %cubeStr=cell2grid(coord,gasStruct.EnergyDissipation(gasMask),cellSize,...
            %    'ngrid',Ngrid,'extensive','weights',mass,'box',boxSize);
            cube=cubeStr.cube.*illUnits.EnergyDissipationUnit; %in K
            logFlag=true;
            %weight=cubeStr.weights; % cube of mass in each uniform grid cell
            weight=ones(size(cube));
            bartag='$\log E_\mathrm{dis}\,[10^{45}\mathrm{erg/yr}]$';
            slTypeDef='avg';
            printTypeTag='ediss';
            
        case {'pressure','press','p'}
            if ~isfield(gasStruct,'Pressure')
                gasStruct=illustris.utils.addPressure(gasStruct); %in K
            end
            
            pre=gasStruct.Pressure(gasMask);
            msk=sfrMask;
            cubeStr=cell2grid(coord(:,msk),pre(msk),cellSize(msk),...
                'ngrid',Ngrid,'intensive','weights',mass(msk),'box',boxSize);
            cube=cubeStr.cube./1e-10; %in 10^-10 erg/cm^3
            logFlag=true;
            weight=cubeStr.weights; % cube of mass in each uniform grid cell
            bartag='$\log P\,[\mathrm{10^{-10}erg/cm^3}]$';
            slTypeDef='avg';
            printTypeTag='press';
            cmap=brewermap(256,'PuRd');
            
        case {'xray','lx'}
        case {'metalicity','metals','ztot'}
            cubeStr=cell2grid(coord,gasStruct.GFM_Metallicity(gasMask),cellSize,...
                'ngrid',Ngrid,'intensive','weights',mass,'box',boxSize);
            cube=cubeStr.cube; %i
            logFlag=true;
            weight=cubeStr.weights; % cube of mass in each uniform grid cell
            bartag='$\log Z\,[\mathrm{Z_\odot}]$';
            slTypeDef='avg';
            cmap=brewermap(256,'BuPu');
        case {'potential','potent'}
        case {'kinetic'}
        case {'velocity','vel'}
            vv=zeros(3,cellCount);
            for k=1:3
                vv(k,:)=gasStruct.Velocities(k,gasMask)-vcm(k);
            end
            vv=sqrt(sum(vv.^2,1));
            cubeStr=cell2grid(coord,vv,cellSize,...
                'ngrid',Ngrid,'intensive','weights',mass,'box',boxSize);
            
            cube=cubeStr.cube.*illustris.utils.velocityFactor(illUnits.snap,'gas'); %in km/sec
            logFlag=true;
            weight=cubeStr.weights; % cube of mass in each uniform grid cell
            bartag='$\log |v|\,[\mathrm{km/sec}]$';
            slTypeDef='avg';
            printTypeTag='vMag';
            
        case {'vr','radialvel'}
            vv=zeros(3,cellCount);
            dist=sqrt(sum(coord.^2,1));
            vr=zeros(3,cellCount);
            for k=1:3
                vv(k,:)=gasStruct.Velocities(k,gasMask)-vcm(k);
                vr(k,:)=vv(k,:).*coord(k,:)./dist;
                
            end
            vrr=sum(vr,1);
            
            cubeStr=cell2grid(coord,vrr,cellSize,...
                'ngrid',Ngrid,'intensive','weights',mass,'box',boxSize);
            cube=cubeStr.cube.*illustris.utils.velocityFactor(illUnits.aexp,'gas'); %in km/sec
            logFlag=false;
            weight=cubeStr.weights; % cube of mass in each uniform grid cell
            bartag='$v_r\,[\mathrm{km/sec}]$';
            slTypeDef='avg';
            printTypeTag='vRad';
            map = brewermap(256,'*RdBu');
            
        otherwise
            error('%s - unknown data type: %s',current_function().upper,type)
    end
end
cube=cube.*normfactor;

boxSize=cubeStr.boxSide;

%% prepare velocity stuff
if velFlag || strmFlag
    
    for k=1:3
        if ~vcubeFlag(k)
            
            
            vv=gasStruct.Velocities(k,gasMask)-vcm(k);
            
            
            vcubeStr=cell2grid(coord,vv,cellSize,...
                'ngrid',Ngrid,'intensive','weights',mass,'box',boxSize);
            switch k
                case 1
                    vx=vcubeStr.cube;
                case 2
                    vy=vcubeStr.cube;
                case 3
                    vz=vcubeStr.cube;
            end
        end
    end
    vWeights=vcubeStr.weights; % cube of mass in each uniform grid cell
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
%             contrWeight=ones(size(contrCube));
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
if velFlag || strmFlag
    [v_x, v_y]=mk_vfield(vx,vy,vz,vWeights,slind);
    clear vx vy vz;
    diluted_len = length(1:dilute:Ngrid);
    diluted_jump = 2.*halfSide/(diluted_len-1);
    notdiluted_jump = 2.*halfSide/(Ngrid-1);
    [xxv, yyv] = meshgrid(-halfSide:diluted_jump:halfSide, -halfSide:diluted_jump:halfSide);
    [xxs, yys] = meshgrid(-halfSide:notdiluted_jump:halfSide,-halfSide:notdiluted_jump:halfSide);
end



%% draw map
if ~clim_Flag
    clims(1) = min(slice(slice>-Inf));
    clims(2) = max(slice(slice<Inf));
    
end
for projection = 1:3
    if plotproj(projection)
        
        if newFigFlag
            if ~exist('hf')
                hf=figure('color',figureColorBk);
            else
                hf(end+1)=figure('color',figureColorBk);
            end
        elseif exist('hf')
            figure(hf)
        else
            error('%s -  No valid Figure handle given',current_function().upper);
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
        
        switch(lower(nanVal))
            case 'none'
            case{'min'}
                msk=isnan(imj);
                imj(msk)=clims(1)-0.1.*abs(clims(1));
            
                case{'max'}
                msk=isnan(imj);
                imj(msk)=clims(2)+0.1.*abs(clims(2));
            otherwise
               error('%s - illegal value for nanVal: %s',current_function().upper,nanVal)
        end
                
        
        %if ~pointPlotFlag
        imagesc(side,side,imj,clims);%axis equal tight;
        
        %         else
        %             cl=diff(side)/Ngrid;
        %             ss=side + 0.5*cl.*[1 -1];
        %             xx=linspace(ss(1),ss(2),Ngrid);
        %
        %             scatter(xx,xx,3,imj)
        
        hold on
        
        
        
        %% add velocity field
        if velFlag || strmFlag
            resampled_v_x = imresize(squeeze(v_x(:,:,projection)), [diluted_len diluted_len], 'bilinear');
            resampled_v_y = imresize(squeeze(v_y(:,:,projection)), [diluted_len diluted_len], 'bilinear');
            if strmFlag
                h = streamslice(xxs,yys,squeeze(v_x(sideInd,sideInd,projection)),squeeze(v_y(sideInd,sideInd,projection)),streamDense);
                set(h,'color',streamColor)
            end
            if velFlag
                quiver(xxv,yyv,resampled_v_x,resampled_v_y, 1.3,'color',vfieldColor);
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
                    circStruct(indC).color=[0 0 0];
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
        
        %% draw a line
        if exist('linStruct','var')
            for indL=1:length(linStruct) % go over all lines
                
                if ~isfield(linStruct(indL),'color')||isempty(linStruct(indL).color) % set color  if none is given
                    linStruct(indL).color=[0 0 0];
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
                        error('%s - Illegal direction in draw line: %s',current_function().upper,linStruct.dir)
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
                    error('%s - Arrow structure must have a START field',current_function().upper)
                end
                if ~isfield(arrowStruct(indA),'stop') % set line type  if none is given
                    error('%s - Arrow structure must have a STOP field',current_function().upper)
                end
                if ~isfield(arrowStruct(indA),'color') || isempty(arrowStruct(indA).color) % set color  if none is given
                    arrowStruct(indA).color=[0 0 0];
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
        
        
        if(comoveFlag && zred>0)
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
        
        
        set(gca,'Ydir','normal','Fontsize',14,...
            'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1],...
            'XLim',side,'YLim',side,'color',figureColorText,...
            'xcolor',figureColorText,'ycolor',figureColorText);    %'XTickLabel',xtickLabels,'YTickLabel',ytickLabels
        % 'XTick',xticks ,'YTick',yticks,... %'TickLength',[-0.015 -0.015],...
        
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
                error('%s - Illegal zoom box',current_function().upper)
            end
            set(gca,'XLim',zb(1,:),'YLim',zb(2,:))
            
            
        end
        
        %% do bar and colormap stuff
        caxis(clims);
        
        if brewerFlag
            map = brewermap(256,brewMap);
        end
        
        if backgroundFlag
            map(1,:)=figureColorBk;
        end
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
        
        set(gcf,'Colormap',map);
        %set(gcf,'Colormap',avijet);
        %title(sprintf('%s %s, Thickness=%s Mpc/h',CLUSTER,type,num2str(thick,3)),'Fontsize',12,'Interpreter','latex');
        
        %% set title
        if titleFlag
            ht=titlemine(titletag);
            set(ht,'color',figureColorText)
        end
        %                     switch titletype
        %                         case 'full'
        %                             mvstr=mk_mvir_string(MVIR);
        %                             titlemine(sprintf('%s %s, %s $z=%3.2g$',CLUSTER,mvstr,titletag,zred));
        %                         case 'custom'
        %                             titlemine(titletag);
        %                         otherwise
        %                             error('MKMAP: Illegal value for title type (full/custom)')
        %                     end
        %                 end
        
        %% add text 
        if exist('textBoxStruct','var')
            for indA=1:length(textBoxStruct)
                if ~isfield(textBoxStruct(indA),'fontsize') || isempty(textBoxStruct(indA).fontsize)% set linewidth  if none is given
                    textBoxStruct(indA).fontsize=20;
                end
                
                if ~isfield(textBoxStruct(indA),'Interpreter') || isempty(textBoxStruct(indA).Interpreter)% set linewidth  if none is given
                    textBoxStruct(indA).Interpreter='latex';
                end
                
                 textHandle=text();
                 setPropertiesStructure(textHandle,textBoxStruct(indA))
                
%                 %             if ~isfield(textBoxStruct(indA),'box') || isempty(textBoxStruct(indA).fontsize)% set linewidth  if none is given
%                 %                     textBoxStruct(indA).fontsize=20;
%                 %             end
%                 
%                 % calculate position of text (given in data units) to position
%                 % within the axes and figure, as needed by annotation (0-1 in
%                 % figure)
%                 axPos=get(gca,'pos');
%                 xl=xlim;
%                 yl=ylim;
%                 pos=textBoxStruct.pos;
%                 
%                 pos(1)=axPos(1)+axPos(3)*(pos(1)-xl(1))./diff(xl);
%                 pos(2)=axPos(2)+axPos(4)*(pos(2)-yl(1))./diff(yl);
%                 pos(3)=axPos(3)*pos(3)./diff(xl);
%                 pos(4)=axPos(4)*pos(4)./diff(yl);
%                 
%                 
%                 
%                 
%                 %textBoxStruct.pos=
%                 annotation(gcf,'textbox',pos,...
%                     'Color',textBoxStruct.color,...
%                     'String', textBoxStruct.str,...
%                     'Interpreter','latex',...
%                     'Fontsize',textBoxStruct.fontsize,...
%                     'FitBoxToText','off','linestyle','none');
                
            end
        end
        
        
        
        
        
        %% add grid
        if gridFlag
            grid;
        end
        
        %% print
        if printFlag
            
            name=sprintf('gasMap%s_%s_%s_snp%i_%s',prjtag,printTypeTag,printtag,illUnits.snap,simDisplayName);
            %name=sprintf('%s/%s_map%s_b%d_%s.%s',printoutdir,CLUSTER,prjtag,boxx,printtag,'%s');
            
            
            printout_fig(gcf,name,'dir',printoutdir,'v',pdfFlag)
            
            
        end
        
        
        
        
    end
end

if saveFigFlag
    
    name=sprintf('%s/figFiles/gasMap_%s_%s_snp%i_%s.fig',DEFAULT_PRINTOUT_DIR,printTypeTag,printtag,illUnits.snap,simDisplayName);
    %name=sprintf('%s/%s_map%s_b%d_%s.%s',printoutdir,CLUSTER,prjtag,boxx,printtag,'%s');
    
    savefig(hf,name,'compact')
    
end

%
%

end

