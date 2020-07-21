function   mkmap2Things(varargin )
%MKMAP plotting the TNG objects
%   Plotting TNG data maps with 2 parameters setting the color/brightness


%% defuals and globals
global DEFAULT_PRINTOUT_DIR
global illUnits
units;

zred=0;
projectiontags={'YZ' 'XZ' 'XY'};
%bartagOveride='';
printoutdir=DEFAULT_PRINTOUT_DIR;
printtag='badtag';
Ngrid=256;
plotproj=[0 0 0]; % plot all three projections
%dilute=4;
%streamDense=1.5;
marks=[];
vcm=[0,0,0];
type='';
thick=-Ngrid   ;
normfactor=[1 1];  % option to enter a normalization for the data
titletag='';
labeltype='full';
titletype='full';
streamColor='k';
vfieldColor='w';
markColor='k';
stlw=1;
slTypeDef{1}='avg'; % is slice an average of slices or sum of them
slTypeDef{2}='avg'; % is slice an average of slices or sum of them

%hf=0;
axesHandle=0;
boxSize=0;


%Flags
cubeFlag=[false false];
typeFlag=[false false];
wtFlag=[false false];
%labelFlag=false;
titleFlag=true;
printFlag=false;
parLim_Flag=[false false]; % no colorbar limits given, use min max
velFlag=false;
strmFlag=false;
%mvirboxFlag=false;
logFlag=[false false];
comoveFlag=true;
vcubeFlag(1:3)=false;
rogFlag=false;
zoomFlag=false;
gridFlag=false;
brewerFlag=false;
%gasStructFlag=false;
newFigFlag=true;
backgroundFlag=false;
cmap = brewermap(256,'*Spectral');
saveFigFlag=false;
colorboxFlag=true;
%% reading arguments

i=1;
while i<=length(varargin)
    switch lower(varargin{i})
        
        % arguments for calculating cube based on gas Structure
        case{'gas','gasstruct'}
            gasStructFlag=true;
            i=i+1;
            gasStruct=varargin{i};
        case{'star','starstruct','stars'}
            %starStructFlag=true;
            i=i+1;
            starStruct=varargin{i};
        case{'main','maintype'} %prepare data cube of this type for main parameter
            i=i+1;
            typeFlag(1)=true;
            typeMain=varargin{i};
        case{'sec','sectype','secondary'} %prepare data cube of this type for sec. parameter
            i=i+1;
            typeFlag(2)=true;
            typeSec=varargin{i};
        case{'ngrid','ng'}  % set the grid size for the uniform cube
            i=i+1;
            Ngrid=varargin{i};
        case{'mask','gasmask'}
            i=i+1;
            gasMask=varargin{i};
        case{'starmask'}
            i=i+1;
            starMask=varargin{i};
        case{'box','boxsize'}
            i=i+1;
            boxSize=varargin{i};
            
            
            
            % arguments for plotting a given cube (main parameter)
        case{'datamain', 'maindata', 'maindatacube','maincube'} %plot a given cube;
            i=i+1;
            cubeFlag(1)=true;
            cubeStr=varargin{i};
            cubeMain=cubeStr.cube;
            printTypeTag='data';
        case{'mainweightcube','mainwtcube','mainwtdatacube','mainwt'} % weight for cube
            i=i+1;
            wtFlag(1)=true;
            weightMain=varargin{i};
            
            % set whether to average or sum up slice
        case{'mainavg','mainaverage','avg1','1avg'}
            slType(1)='avg';
        case{'mainsum','sum1','1sum'}
            slType(1)='sum';
            
            % arguments for plotting a given cube (secondary parameter)
        case{'datasec', 'secata', 'secdatacube','seccube'} %plot a given cube;
            i=i+1;
            cubeFlag(2)=true;
            cubeStr=varargin{i};
            cubeSec=cubeStr.cube;
            %printTypeTag='data';
        case{'secweightcube','secwtcube','secwtdatacube','secwt'} % weight for cube
            i=i+1;
            wtFlag(2)=true;
            weightSec=varargin{i};
            
            % set whether to average or sum up slice
        case{'secavg','secaverage','avg2','2avg'}
            slType(2)='avg';
        case{'secsum','sum2','2sum'}
            slType(2)='sum';
            
            
            
            
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
                error('MKMAPGAS - vcm must be a 3-component array')
            end
            
            % basic arguments for plotting
        case{'logmain','mainlog'}
            logFlag(1)=true;
        case{'nologmain','mainlinear','linearmain'}
            logFlag(1)=false;
            
        case{'logsec','seclog'}
            logFlag(2)=true;
        case{'nologsec','seclinear','linearsec'}
            logFlag(2)=false;
            
            
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
        case {'mainlimits','mainlim','limmain','limitsmain'}
            i=i+1;
            parLim_Flag(1)=true;
            limMain=varargin{i};
            if diff(limMain)==0
                parLim_Flag(1)=false;
            end
        case {'seclimits','seclim','limsec','limitssec'}
            i=i+1;
            parLim_Flag(2)=true;
            limSec=varargin{i};
            if diff(limSec)==0
                parLim_Flag(2)=false;
            end
            
            
        case {'normalizemain','normmain','factormain','normfactormain', 'mainnorm','norm1','mainfactor','factor1','mainnormfactor'}
            i=i+1;
            normfactor(1)=varargin{i};
        case {'normalizesec','normsec','factorsec','normfactorsec', 'secnorm','norm2','secfactor','factor2','secnormfactor'}
            i=i+1;
            normfactor(2)=varargin{i};
            
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
        case {'brewer','brewermap'}
            i=i+1;
            brewMap=varargin{i};
            brewerFlag=true;
        case {'cmap','colormap'}
            i=i+1;
            cmap=varargin{i};
            brewerFlag=false;
            if size(cmap,2)~=3
                error('mkmapGas: Illegal colormap')
            end
            
            
        case {'nocolorbox','nocolorbar','nobar'}
                colorboxFlag=false;
        case {'nobackground'}
            backgroundFlag=false;
%         case 'xticks'
%             i=i+1;
%             xticks=varargin{i};
%         case 'yticks'
%             i=i+1;
%             yticks=varargin{i};
%         case 'fullticks'
%             i=i+1;
%             tickStruct=varargin{i};
%             if ~isstruct(tickStruct)
%                 error('MKMAP - fullTick input must be a structure')
%             end
        case {'labels'}
            i=i+1;
            labeltype=varargin{i};
            %labelFlag=(strcmp(l,'yes') || strcmp(answer,'on') || strcmp(answer));
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
            if length(bartagOveride)~=2
                error('MKMAP2THinggs - there shoulwd be only 2 colorbar tags')
            end
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
            
            
            %% printing
        case {'print'}
            printFlag=true;
            i=i+1;
            printtag=varargin{i};
        case {'outputdir','printout','printoutdir'}
            i=i+1;
            printoutdir=varargin{i};
        case {'savefig'}
            i=i+1;
            printtag=varargin{i};
            saveFigFlag=true;
        otherwise
            error('mkmapGas: illegal argument %s',varargin{i})
    end
    i=i+1;
end


if ~any(plotproj)
    plotproj=[1 1 1];
end

% if ~cubeFlag && ~(gasStructFlag && typeFlag)
%     error('mkmapGas: must enter both gas struct and data type to plot');
% end


% if ~(cubeFlag || typeFlag)
%     error('mkmapGas: must enter data or datatype')
% end
% if (cubeFlag && typeFlag)
%     error('mkmapGas: too many data arguments')
% end


if ~wtFlag(1)
    weightMain=ones(Ngrid,Ngrid,Ngrid);
end

if ~wtFlag(2)
    weightSec=ones(Ngrid,Ngrid,Ngrid);
end


%% add missing fields
if (gasStructFlag)
    
    if ~isfield(gasStruct,'newCoord')
        error('MKMAPGAS - object has not been centered !!')
    end
end


%% build Main parameter data cube
if typeFlag(1)
    switch lower(typeMain)
        case {'flux','massflux','f'}
            %             %fluxnorm=(0.1117.*(MVIR./1e15).^0.15*(1+zred)^2.25)./(4.*pi); % normalization for flux
            %
            %             % find vr
            %             vv=zeros(3,gasStruct.count);
            %             dist=sqrt(sum(gasStruct.newCoord.^2,1));
            %             vr=zeros(3,gasStruct.count);
            %             for k=1:3
            %                 vv(k,:)=gasStruct.Velocities(k,:)-vcm(k);
            %                 vr(k,:)=vv(k,:).*gasStruct.newCoord(k,:)./dist;
            %
            %             end
            %             vrr=sum(vr,1);
            %
            %             cubeStr=cell2grid(gasStruct.newCoord,vrr.*gasStruct.Masses,(gasStruct.Masses./gasStruct.Density).^(1/3),...
            %                 'ngrid',Ngrid,'extensive','box',boxSize);
            %             cube=cubeStr.cube./(cubeStr.cellVol).*illUnits.densityUnit;  % Msun/kpc^3*km/sec
            %             logFlag=false;
            %             %weight=cubeStr.weights; % cube of mass in each uniform grid cell
            %             bartag='$v_r\,[\mathrm{km/sec}]$';
            %             slTypeDef='sum';
            %
            %
            %             %             %find radial velocity component
            %             %             [meshY, meshX, meshZ] = meshgrid(1:size(vy,1), 1:size(vx,2), 1:size(vz,3));
            %             %             %convert to center origin coordinates
            %             %             meshX = meshX - (size(vx,1)+1)/2 -cm(1);
            %             %             meshY = meshY - (size(vy,2)+1)/2 -cm(2);
            %             %             meshZ = meshZ - (size(vz,3)+1)/2 -cm(3);
            %             %             % Fix Units (to be in Mpc)
            %             %             meshX = meshX * ((boxx/hub)/NCELL);
            %             %             meshY = meshY * ((boxx/hub)/NCELL);
            %             %             meshZ = meshZ * ((boxx/hub)/NCELL);
            %             %
            %             %             rcube2=meshX.^2+meshY.^2+meshZ.^2 ; % r^2 cube in Mpc
            %             %             vr=(vx.*meshX+vy.*meshY+vz.*meshZ)./sqrt(rcube2) ; %radial velocity
            %             %
            %             %             switch SIMTYPE
            %             %                 case 'CSF'
            %             %                     [Mg200 , ~, ~]=read_Mass_Profiles(RVIR);
            %             %                 case 'Adiabatic'
            %             %                     bx=2^(ceil(log(2*get_rvir*hub)/log(2)));
            %             %                     mg=RHOG(bx).*(bx./NCELL./hub).^3;
            %             %                     mgP=MAKE_CUM_PROFILE_FROM_CUBE(mg);
            %             %                     ind=ceil(get_rvir/(0.5*bx/hub/NCELL));
            %             %                     Mg200=mgP(ind);
            %             %             end
            %             %
            %             %             cube=rog.*vr.*rcube2./Mg200.*(km/Mpc*Gyr)./fluxnorm ;
            %             %             weight=ones(size(cube));
            %             %             bartag='$\frac{\dot{M}_{gas}}{d\Omega}\,[\frac{\dot{M}_{gas}}{d\Omega}|_\mathrm{vir} ]$';
            %             %             %bartag='$\frac{\dot{M}/M_{gas}}{d\Omega} [\frac{\dot M}{M_{gas}}\left|_{\mathrm{virial}}]$';
            %             %             %bartag='$\dot{M}_{gas}$';
            %             %
            %             %             printTypeTag='flux';
            %             %
            %             %             clear vr rcube2 meshX meshY meshZ
            %
            %
        case {'density','gasdensity','rho','gas','rhogas'}
            cubeStr=cell2grid(gasStruct.newCoord,gasStruct.Masses,(gasStruct.Masses./gasStruct.Density).^(1/3),'ngrid',Ngrid,'extensive','box',boxSize);
            
            %cube=(cubeStr.cube.*massUnit)./(cubeStr.cellVol.*lengthUnit^3);
            cubeMain=(cubeStr.cube)./(cubeStr.cellVol).*illUnits.densityUnit; % in Msun/kpc^3
            weightMain=ones(size(cubeMain));
            logFlag(1)=true;
            bartag{1}='$\log \rho_\mathrm{g}\,[\mathrm{M_\odot/kpc^3}]$';
            slTypeDef{1}='avg';
            printTypeTag{1}='dens';
            
        case {'numberdensity','n'}
            cubeStr=cell2grid(gasStruct.newCoord,gasStruct.Masses,(gasStruct.Masses./gasStruct.Density).^(1/3),'ngrid',Ngrid,'extensive','box',boxSize);
            
            %cube=(cubeStr.cube.*massUnit)./(cubeStr.cellVol.*lengthUnit^3);
            
            %rhoFac=illUnits.densityUnit.*(Units.Ms/Units.kpc^3/Units.mm); %in cm^-3
            cubeMain=(cubeStr.cube)./(cubeStr.cellVol).*illUnits.numberDensityFactor; % in cm^-3
            weightMain=ones(size(cubeMain));
            logFlag(1)=true;
            bartag{1}='$\log n\,[\mathrm{cm^{-3}}]$';
            slTypeDef{1}='avg';
            
            printTypeTag{1}='nDens';
            
            
        case {'entropy','ent','k','s'}
            if ~isfield(gasStruct,'Entropy')
                gasStruct=illustris.utils.addEntropy(gasStruct); %in KeV cm^2
            end
            
            cubeStr=cell2grid(gasStruct.newCoord,gasStruct.Entropy,(gasStruct.Masses./gasStruct.Density).^(1/3),...
                'ngrid',Ngrid,'intensive','weights',gasStruct.Masses,'box',boxSize);
            cubeMain=cubeStr.cube; %in K
            logFlag(1)=true;
            weightMain=cubeStr.weights; % cube of mass in each uniform grid cell
            bartag{1}='$\log S\,[\mathrm{KeV\, cm^2}]$';
            slTypeDef{1}='avg';
            printTypeTag{1}='ent';
            
        case {'temperature','temp','t','tmp'}
            
            if ~isfield(gasStruct,'Temperature')
                gasStruct.Temperature=illustris.utils.calcTemperature(gasStruct.InternalEnergy,gasStruct.ElectronAbundance); %in K
            end
            
            cubeStr=cell2grid(gasStruct.newCoord,gasStruct.Temperature,(gasStruct.Masses./gasStruct.Density).^(1/3),...
                'ngrid',Ngrid,'intensive','weights',gasStruct.Masses,'box',boxSize);
            cubeMain=cubeStr.cube; %in K
            logFlag(1)=true;
            weightMain=cubeStr.weights; % cube of mass in each uniform grid cell
            bartag{1}='$\log T\,[\mathrm{K}]$';
            slTypeDef{1}='avg';
            printTypeTag{1}='temp';
        case {'mach'}
        case {'shock','shockenergy','energydissipation','dissipation'}
            
        case {'pressure','press','p'}
        case {'xray','lx'}
        case {'metalicity','metals','ztot'}
        case {'potential','potent'}
        case {'kinetic'}
        case {'velocity','vel'}
            vv=zeros(3,gasStruct.count);
            for k=1:3
                vv(k,:)=gasStruct.Velocities(k,:)-vcm(k);
            end
            vv=sqrt(sum(vv.^2,1));
            cubeStr=cell2grid(gasStruct.newCoord,vv,(gasStruct.Masses./gasStruct.Density).^(1/3),...
                'ngrid',Ngrid,'intensive','weights',gasStruct.Masses,'box',boxSize);
            cubeMain=cubeStr.cube; %in K
            logFlag(1)=true;
            weightMain=cubeStr.weights; % cube of mass in each uniform grid cell
            bartag{1}='$\log |v|\,[\mathrm{km/sec}]$';
            slTypeDef{1}='avg';
            printTypeTag{1}='vMag';
            
        case {'vr','radialvel'}
            vv=zeros(3,gasStruct.count);
            dist=sqrt(sum(gasStruct.newCoord.^2,1));
            vr=zeros(3,gasStruct.count);
            for k=1:3
                vv(k,:)=gasStruct.Velocities(k,:)-vcm(k);
                vr(k,:)=vv(k,:).*gasStruct.newCoord(k,:)./dist;
                
            end
            vrr=sum(vr,1);
            
            cubeStr=cell2grid(gasStruct.newCoord,vrr,(gasStruct.Masses./gasStruct.Density).^(1/3),...
                'ngrid',Ngrid,'intensive','weights',gasStruct.Masses,'box',boxSize);
            cubeMain=cubeStr.cube; %in K
            logFlag(1)=false;
            weightMain=cubeStr.weights; % cube of mass in each uniform grid cell
            bartag{1}='$v_r\,[\mathrm{km/sec}]$';
            slTypeDef{1}='avg';
            printTypeTag{1}='vRad';
            
            
        otherwise
            error('mkmapGas - unknown main data type: %s',typeMain)
    end
end
cubeMain=cubeMain.*normfactor(1);


%% build Secondary  parameter data cube
if typeFlag(2)
    switch lower(typeSec)
        case {'flux','massflux','f'}
            %             %fluxnorm=(0.1117.*(MVIR./1e15).^0.15*(1+zred)^2.25)./(4.*pi); % normalization for flux
            %
            %             % find vr
            %             vv=zeros(3,gasStruct.count);
            %             dist=sqrt(sum(gasStruct.newCoord.^2,1));
            %             vr=zeros(3,gasStruct.count);
            %             for k=1:3
            %                 vv(k,:)=gasStruct.Velocities(k,:)-vcm(k);
            %                 vr(k,:)=vv(k,:).*gasStruct.newCoord(k,:)./dist;
            %
            %             end
            %             vrr=sum(vr,1);
            %
            %             cubeStr=cell2grid(gasStruct.newCoord,vrr.*gasStruct.Masses,(gasStruct.Masses./gasStruct.Density).^(1/3),...
            %                 'ngrid',Ngrid,'extensive','box',boxSize);
            %             cube=cubeStr.cube./(cubeStr.cellVol).*illUnits.densityUnit;  % Msun/kpc^3*km/sec
            %             logFlag=false;
            %             %weight=cubeStr.weights; % cube of mass in each uniform grid cell
            %             bartag='$v_r\,[\mathrm{km/sec}]$';
            %             slTypeDef='sum';
            %
            %
            %             %             %find radial velocity component
            %             %             [meshY, meshX, meshZ] = meshgrid(1:size(vy,1), 1:size(vx,2), 1:size(vz,3));
            %             %             %convert to center origin coordinates
            %             %             meshX = meshX - (size(vx,1)+1)/2 -cm(1);
            %             %             meshY = meshY - (size(vy,2)+1)/2 -cm(2);
            %             %             meshZ = meshZ - (size(vz,3)+1)/2 -cm(3);
            %             %             % Fix Units (to be in Mpc)
            %             %             meshX = meshX * ((boxx/hub)/NCELL);
            %             %             meshY = meshY * ((boxx/hub)/NCELL);
            %             %             meshZ = meshZ * ((boxx/hub)/NCELL);
            %             %
            %             %             rcube2=meshX.^2+meshY.^2+meshZ.^2 ; % r^2 cube in Mpc
            %             %             vr=(vx.*meshX+vy.*meshY+vz.*meshZ)./sqrt(rcube2) ; %radial velocity
            %             %
            %             %             switch SIMTYPE
            %             %                 case 'CSF'
            %             %                     [Mg200 , ~, ~]=read_Mass_Profiles(RVIR);
            %             %                 case 'Adiabatic'
            %             %                     bx=2^(ceil(log(2*get_rvir*hub)/log(2)));
            %             %                     mg=RHOG(bx).*(bx./NCELL./hub).^3;
            %             %                     mgP=MAKE_CUM_PROFILE_FROM_CUBE(mg);
            %             %                     ind=ceil(get_rvir/(0.5*bx/hub/NCELL));
            %             %                     Mg200=mgP(ind);
            %             %             end
            %             %
            %             %             cube=rog.*vr.*rcube2./Mg200.*(km/Mpc*Gyr)./fluxnorm ;
            %             %             weight=ones(size(cube));
            %             %             bartag='$\frac{\dot{M}_{gas}}{d\Omega}\,[\frac{\dot{M}_{gas}}{d\Omega}|_\mathrm{vir} ]$';
            %             %             %bartag='$\frac{\dot{M}/M_{gas}}{d\Omega} [\frac{\dot M}{M_{gas}}\left|_{\mathrm{virial}}]$';
            %             %             %bartag='$\dot{M}_{gas}$';
            %             %
            %             %             printTypeTag='flux';
            %             %
            %             %             clear vr rcube2 meshX meshY meshZ
            %
            %
        case {'density','gasdensity','rho','gas','rhogas'}
            cubeStr=cell2grid(gasStruct.newCoord,gasStruct.Masses,(gasStruct.Masses./gasStruct.Density).^(1/3),'ngrid',Ngrid,'extensive','box',boxSize);
            
            %cube=(cubeStr.cube.*massUnit)./(cubeStr.cellVol.*lengthUnit^3);
            cubeSec=(cubeStr.cube)./(cubeStr.cellVol).*illUnits.densityUnit; % in Msun/kpc^3
            weightSec=ones(size(cubeSec));
            logFlag(2)=true;
            bartag{2}='$\log \rho_\mathrm{g}\,[\mathrm{M_\odot/kpc^3}]$';
            slTypeDef{2}='avg';
            printTypeTag{2}='dens';
            
        case {'numberdensity','n','number'}
            cubeStr=cell2grid(gasStruct.newCoord,gasStruct.Masses,(gasStruct.Masses./gasStruct.Density).^(1/3),'ngrid',Ngrid,'extensive','box',boxSize);
            
            %cube=(cubeStr.cube.*massUnit)./(cubeStr.cellVol.*lengthUnit^3);
            
            %rhoFac=illUnits.densityUnit.*(Units.Ms/Units.kpc^3/Units.mm); %in cm^-3
            cubeSec=(cubeStr.cube)./(cubeStr.cellVol).*illUnits.numberDensityFactor; % in cm^-3
            weightSec=ones(size(cubeSec));
            logFlag(2)=true;
            bartag{2}='$\log n\,[\mathrm{cm^{-3}}]$';
            slTypeDef{2}='avg';
            
            printTypeTag{2}='nDens';
            
            
        case {'entropy','ent','k','s'}
            if ~isfield(gasStruct,'Entropy')
                gasStruct=illustris.utils.addEntropy(gasStruct); %in KeV cm^2
            end
            
            cubeStr=cell2grid(gasStruct.newCoord,gasStruct.Entropy,(gasStruct.Masses./gasStruct.Density).^(1/3),...
                'ngrid',Ngrid,'intensive','weights',gasStruct.Masses,'box',boxSize);
            cubeSec=cubeStr.cube; %in K
            logFlag(2)=true;
            weightSec=cubeStr.weights; % cube of mass in each uniform grid cell
            bartag{2}='$\log S\,[\mathrm{KeV\, cm^2}]$';
            slTypeDef{2}='avg';
            printTypeTag{2}='ent';
            
        case {'temperature','temp','t','tmp'}
            
            if ~isfield(gasStruct,'Temperature')
                gasStruct.Temperature=illustris.utils.calcTemperature(gasStruct.InternalEnergy,gasStruct.ElectronAbundance); %in K
            end
            
            cubeStr=cell2grid(gasStruct.newCoord,gasStruct.Temperature,(gasStruct.Masses./gasStruct.Density).^(1/3),...
                'ngrid',Ngrid,'intensive','weights',gasStruct.Masses,'box',boxSize);
            cubeSec=cubeStr.cube; %in K
            logFlag(2)=true;
            weightSec=cubeStr.weights; % cube of mass in each uniform grid cell
            bartag{2}='$\log T\,[\mathrm{K}]$';
            slTypeDef{2}='avg';
            printTypeTag{2}='temp';
        case {'mach'}
        case {'shock','shockenergy','energydissipation','dissipation'}
            
        case {'pressure','press','p'}
        case {'xray','lx'}
        case {'metalicity','metals','ztot'}
        case {'potential','potent'}
        case {'kinetic'}
        case {'velocity','vel'}
            vv=zeros(3,gasStruct.count);
            for k=1:3
                vv(k,:)=gasStruct.Velocities(k,:)-vcm(k);
            end
            vv=sqrt(sum(vv.^2,1));
            cubeStr=cell2grid(gasStruct.newCoord,vv,(gasStruct.Masses./gasStruct.Density).^(1/3),...
                'ngrid',Ngrid,'intensive','weights',gasStruct.Masses,'box',boxSize);
            cubeSec=cubeStr.cube; %in K
            logFlag(2)=true;
            weightSec=cubeStr.weights; % cube of mass in each uniform grid cell
            bartag{2}='$\log |v|\,[\mathrm{km/sec}]$';
            slTypeDef{2}='avg';
            printTypeTag{2}='vMag';
            
        case {'vr','radialvel'}
            vv=zeros(3,gasStruct.count);
            dist=sqrt(sum(gasStruct.newCoord.^2,1));
            vr=zeros(3,gasStruct.count);
            for k=1:3
                vv(k,:)=gasStruct.Velocities(k,:)-vcm(k);
                vr(k,:)=vv(k,:).*gasStruct.newCoord(k,:)./dist;
                
            end
            vrr=sum(vr,1);
            
            cubeStr=cell2grid(gasStruct.newCoord,vrr,(gasStruct.Masses./gasStruct.Density).^(1/3),...
                'ngrid',Ngrid,'intensive','weights',gasStruct.Masses,'box',boxSize);
            cubeSec=cubeStr.cube; %in K
            logFlag(2)=false;
            weightSec=cubeStr.weights; % cube of mass in each uniform grid cell
            bartag{2}='$v_r\,[\mathrm{km/sec}]$';
            slTypeDef{2}='avg';
            printTypeTag{2}='vRad';
            
            
        otherwise
            error('mkmapGas - unknown secondary data type: %s',typeSec)
    end
end
cubeSec=cubeSec.*normfactor(2);


%% prepare velocity stuff
if velFlag || strmFlag
    
    for k=1:3
        if ~vcubeFlag(k)
            
            
            vv=gasStruct.Velocities(k,:)-vcm(k);
            
            
            vcubeStr=cell2grid(gasStruct.newCoord,vv,(gasStruct.Masses./gasStruct.Density).^(1/3),...
                'ngrid',Ngrid,'intensive','weights',gasStruct.Masses,'box',boxSize);
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

%% prepare slice

%define slice index
if thick>0
    thk=ceil(0.5.*thick./cubeStr.boxSide.*cubeStr.Ngrid);   %% thick is in comoving Mpc/h
else
    thk=0.5*abs(thick);   %% defualt value of 6 cells for slice
end

slind=(floor(0.5.*Ngrid)-thk+1):1:(ceil(0.5.*Ngrid)+thk);

sideInd=1:Ngrid;
halfSide=0.5*cubeStr.boxSide;

side=[-halfSide halfSide];

% if(~comoveFlag)
%     side=side./(1+zred);
%     halfSide=halfSide./(1+zred);
% end

if ~exist('slType','var') % choose whether to average or sum slice
    slType=slTypeDef;
end


%% perpare Main slice to be plotted.


sliceMain=mk_slice(cubeMain,weightMain,slind,slType{1});
if logFlag(1)
    sliceMain=log10(sliceMain);
end

sliceSec=mk_slice(cubeSec,weightSec,slind,slType{2});
if logFlag(2)
    sliceSec=log10(sliceSec);
end


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



%% get limits if not defined
if ~parLim_Flag(1)
    limMain(1) = min(sliceMain(sliceMain>-Inf));
    limMain(2) = max(sliceMain(sliceMain<Inf));
end

if ~parLim_Flag(2)
    limSec(1) = min(sliceSec(sliceSec>-Inf));
    limSec(2) = max(sliceSec(sliceSec<Inf));
end

for projection = 1:3
    if plotproj(projection)
        
        if newFigFlag
            if ~exist('hf')
                hf=figure;
            else
                hf(end+1)=figure;
            end
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
        
        %% prepare color map
        if brewerFlag
            cmap = brewermap(256,brewMap);
        end
        if backgroundFlag
            cmap(1,:)=[1 1 1];
        end
        
        
        %% prepare 2-Par map
        res=mkmap2par(sliceMain(:,:,projection),sliceSec(:,:,projection),cmap,limMain,limSec);
        
        
        %% start plotting
        axes1 = axes('Parent',hf(end));
        imagesc(side,side,res.rgbOut); %axis equal tight;
        set(gca,'Ydir','Normal','Fontsize',12,'XLim',side,'YLim',side); 
           % 'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1],...
               %'XTickLabel',xtickLabels,'YTickLabel',ytickLabels
        % 'XTick',xticks ,'YTick',yticks,... %'TickLength',[-0.015 -0.015],...)
        
        
        
        
        
        
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
        
        
        
        
        %% draw circles
        if exist('circStruct','var')
            for indC=1:length(circStruct) % go over all lines
                
                if ~isfield(circStruct(indC),'color')||isempty(circStruct(indC).color) % set color  if none is given
                    circStruct(indC).color='k';
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
                    'EdgeColor',arrowStruct(indA).edgeColor)
            end
            
        end
        hold off
        
        %% put axis labels
        
        prjtag=projectiontags{projection};
        if(comoveFlag && zred>0)
            lunit='[\mathrm{kpc/h\ comoving\,}]';
        else
            lunit='[\mathrm{kpc/h}]';
        end
        
        switch labeltype
            case {'yes','on','full'}
                switch projection
                    case 1  % YZ
                        xlabelmine(sprintf('$Z %s $',lunit));
                        ylabelmine(sprintf('$Y %s $',lunit));
                    case 2  % ZX
                        xlabelmine(sprintf('$X %s $',lunit));
                        ylabelmine(sprintf('$Z %s $',lunit));
                    case 3  % XY
                        xlabelmine(sprintf('$X %s $',lunit));
                        ylabelmine(sprintf('$Y %s $',lunit));
                end
            case {'half','units'}
                xlabelmine(sprintf('$ %s $',lunit));
                ylabelmine(sprintf('$ %s $',lunit));
                
        end
        
      
        
        axes(axes1);
        
        
%         set(gca,'Ydir','normal','Fontsize',14,...
%             'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1],...
%             'XLim',side,'YLim',side);    %'XTickLabel',xtickLabels,'YTickLabel',ytickLabels
        % 'XTick',xticks ,'YTick',yticks,... %'TickLength',[-0.015 -0.015],...
        
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
        
        
        
        %% add grid
        if gridFlag
            grid;
        end
        
          %% add colorbar inset
        if(colorboxFlag)
            axes2 = axes('Parent',hf(end),'Position',[0.12 0.12 0.2 0.2]);
            hold(axes2,'on');
            
            imagesc(res.mainParLim,res.secParLim,res.rgbColorPlane);
            set(gca,'Ydir','Normal','XLim',res.mainParLim,'YLim',res.secParLim)
            xlabel(bartag{1},'Interpreter','latex','Color',[1 1 1],'FontSize',12);
            ylabel(bartag{2},'Interpreter','latex','Color',[1 1 1],'FontSize',12);
            set(axes2,'Layer','top','XAxisLocation','top','XColor',[1 1 1],...
                'YAxisLocation','right','YColor',[1 1 1],'ZColor',[1 1 1]);
            box('on')
        end
        
        
        
        %% print
        if printFlag
            global simDisplayName
            name=sprintf('map%s_%s_%s_%s',prjtag,printTypeTag,printtag,simDisplayName);
            %name=sprintf('%s/%s_map%s_b%d_%s.%s',printoutdir,CLUSTER,prjtag,boxx,printtag,'%s');
            
            
            printout_fig(gcf,name,'dir',printoutdir,'v')
            
            
        end
        
        
        
        
    end
end

if saveFigFlag
    global simDisplayName
    name=sprintf('%s/figFiles/map_%s_%s_%s_%s.fig',printoutdir,printTypeTag{1},printTypeTag{2},printtag,simDisplayName);
    %name=sprintf('%s/%s_map%s_b%d_%s.%s',printoutdir,CLUSTER,prjtag,boxx,printtag,'%s');
    
    savefig(hf,name,'compact')
    
end

%
%

end

