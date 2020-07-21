function mkmap(boxx,varargin)

% basic function for making a map of a cluster plus a velocity field.
% boxx -which boxx to map
% plotproj - which projection to plot, logical vector
% type - which value do we plotve
%clims - limits of colorbar
%velFlag - Flag to create velocity field or not.
% vcm - center of mass velocity should be calculated in advance
% thickf - the thickness of the slice in Mpc/h comoving
% marks - a vector of radii for drawing circles in units of Mpc
% dilute - vector field dilution. a value <0 gives the defualt dilute=4

%% preamble

% argument test
% if nargin < 4
%     error('mkmap: not enough arguments')
%     %elseif nargin > 10
%     %    error('mkmap: too many arguments')
% end

units;

global SIMTYPE
global CLUSTER;
global VCM
global aexpn;
global NCELL;
global zred
global hub;
global DEFAULT_PRINTOUT_DIR
%load(sprintf('%s/virial%d_%s', HALO_PATH, 8,aexpn));
projectiontags={'YZ' 'XZ' 'XY'};
bartagOveride='';

RVIR=get_rvir();
MVIR=get_mvir();
%VVIR=get_vvir();
TVIR=get_tvir();

% defualt values
%boxx=8;
plotproj=[0 0 0]; % plot all three projections
dilute=4;
streamDense=1.5;
marks=[];
vcm=[];
printoutdir=DEFAULT_PRINTOUT_DIR;
printtag='badtag';
printformat='pdf'; %[1 1];  % print eps and png
type='';
thick=-6   ;
normfactor=1;  % option to enter a normalization for the data
titletag='';
labeltype='full';
titletype='full';
rvin=logical([0 1 0]); % must be 3 component vector of 0 or 1 for plotting r500 rvir r200.
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
contrWeight=ones(NCELL,NCELL,NCELL);
alfa=1.0; % parameter for rps plotting.
% default is to plot rvir alone.

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
gridFlag=false;
brewerFlag=false;

%load('MyColormaps','newJet');
%load('MyColormaps','avijet');
map = brewermap(256,'*RdBu');
%map=newJet;
%% reading arguments

i=1;
while i<=length(varargin)
    switch lower(varargin{i})
        case {'boxx','box'}   %which box to plot
            %           boxx=varargin{i};
        case{'data','datacube'} %plot a given cube;
            i=i+1;
            cubeFlag=true;
            cube=varargin{i};
            printTypeTag='data';
        case{'weightcube','wtcube','wtdatacube','wt'} % weight for cube
            i=i+1;
            wtFlag=true;
            weight=varargin{i};
        case{'avg','average'}
            slType='avg';
        case{'sum'}
            slType='sum';
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
        case{'datatype','type'} %prepare data cube of this type
            i=i+1;
            typeFlag=true;
            type=varargin{i};
        case{'log','log10'}
            logFlag=true;
        case{'nolog','linear'}
            logFlag=false;
        case{'width','thick'} %thicknes of slice in Mpc/h comoving (like box)
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
        case {'alfa','alpha'}
            i=i+1;
            alfa=varargin{i};
        case {'normalize','norm'}
            i=i+1;
            normfactor=varargin{i};
        case {'marks','circles'}
            i=i+1;
            marks=varargin{i};
        case 'markcolor'
            i=i+1;
            markColor=varargin{i};
        case {'rv','rv_circles'}   % must be 3 component vector of 0 or 1 for plotting r500 rvir r200.
            i=i+1;
            rvin=logical(varargin{i});
        case {'grid'}
            gridFlag=true;
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
        case 'circle'
             i=i+1;
             circleStruct=varargin{i};
             if ~isstruct(circleStruct)
                 error('MKMAP - circle input must be a structure')
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
                error('mkmap: Illegal colormap')
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
        case {'print'}
            printFlag=true;
        case {'noprint'}
            printFlag=false;
        case {'outputdir','printout','printoutdir'}
            i=i+1;
            printoutdir=varargin{i};
        case {'png'}
            printformat='png';
%         case {'format','printformat'}
%             i=i+1;
%             answer=varargin{i};
%             switch answer
%                 case {'png','PNG'}
%                     printformat=[1 0];
%                 case {'EPS','eps'}
%                     printformat=[0 1];
%                 case {'both'}
%                     printformat=[1 1];
%                 otherwise
%                     error('mkmap printformat: illegal argument %s',varargin{i})
%             end
            
        case {'printtag','printag'}
            i=i+1;
            printtag=varargin{i};
        otherwise
            error('mkmap: illegal argument %s',varargin{i})
    end
    i=i+1;
end

if ~any(plotproj)
    plotproj=[1 1 1];
end

if ~any(boxx==[1 2 4 8])
    error('mkmap: illegal box')
end

if ~(cubeFlag || typeFlag)
    error('mkmap: must enter data or datatype')
end
if (cubeFlag && typeFlag)
    error('mkmap: too many data arguments')
end

if contrFlag
    if ( isempty(contrCube) && isempty(contrType));
        error('mkmap: must enter contoure data or type')
    end
    
    if (~isempty(contrCube) && ~isempty(contrType));
        error('mkmap: too many contour data arguments')
    end
end

if ~rogFlag
    rog=RHOG(boxx);
    rogFlag=true;
end

if (length(vcm)~=3)
    vcm = VCM;
end

if ~wtFlag
    weight=ones(NCELL,NCELL,NCELL);
end

%bcson=sqrt(1.52e8.*T(boxx))./1e5; %speed of sound in km/sec
velFlag2=velFlag || strmFlag ;

if(velFlag2 || strcmpi(type,'flux') || strcmpi(type,'f'))
    if ~rogFlag
        rog=RHOG(boxx);
        rogFlag=true;
    end
    
    [hubX,hubY,hubZ] = hubble_flow(boxx,[0,0,0]);
    vx = Vx(boxx)+hubX-vcm(1);
    vy = Vy(boxx)+hubY-vcm(2);
    vz = Vz(boxx)+hubZ-vcm(3);
    vcubeFlag(:)=true;
end




cm=[0,0,0];
if typeFlag
    switch lower(type)
        case {'flux','massflux','f'}
            fluxnorm=(0.1117.*(MVIR./1e15).^0.15*(1+zred)^2.25)./(4.*pi); % normalization for flux
            
            %find radial velocity component
            [meshY, meshX, meshZ] = meshgrid(1:size(vy,1), 1:size(vx,2), 1:size(vz,3));
            %convert to center origin coordinates
            meshX = meshX - (size(vx,1)+1)/2 -cm(1);
            meshY = meshY - (size(vy,2)+1)/2 -cm(2);
            meshZ = meshZ - (size(vz,3)+1)/2 -cm(3);
            % Fix Units (to be in Mpc)
            meshX = meshX * ((boxx/hub)/NCELL);
            meshY = meshY * ((boxx/hub)/NCELL);
            meshZ = meshZ * ((boxx/hub)/NCELL);
            
            rcube2=meshX.^2+meshY.^2+meshZ.^2 ; % r^2 cube in Mpc
            vr=(vx.*meshX+vy.*meshY+vz.*meshZ)./sqrt(rcube2) ; %radial velocity
            
            switch SIMTYPE
                case 'CSF'
                    [Mg200 , ~, ~]=read_Mass_Profiles(RVIR);
                case 'Adiabatic'
                    bx=2^(ceil(log(2*get_rvir*hub)/log(2)));
                    mg=RHOG(bx).*(bx./NCELL./hub).^3;
                    mgP=MAKE_CUM_PROFILE_FROM_CUBE(mg);
                    ind=ceil(get_rvir/(0.5*bx/hub/NCELL));
                    Mg200=mgP(ind);
            end
            
            cube=rog.*vr.*rcube2./Mg200.*(Units.km/Units.Mpc*Units.Gyr)./fluxnorm ;
            weight=ones(size(cube));
            bartag='$\frac{\dot{M}_{gas}}{d\Omega}\,[\frac{\dot{M}_{gas}}{d\Omega}|_\mathrm{vir} ]$';
            %bartag='$\frac{\dot{M}/M_{gas}}{d\Omega} [\frac{\dot M}{M_{gas}}\left|_{\mathrm{virial}}]$';
            %bartag='$\dot{M}_{gas}$';
            
            printTypeTag='flux';
            
            clear vr rcube2 meshX meshY meshZ
            
        case {'density','gasdensity','rho','gas','rhogas'}
            cube=rog;
            weight=ones(size(cube));
            logFlag=true; %cube=log10(cube);
            bartag='$\log \rho_{gas}\,[\mathrm{M_\odot/Mpc^3}]$';
            printTypeTag='density';
        case {'numberdensity','n'}
            fac=Units.Ms/Units.mm/Units.Mpc^3;
            cube=rog.*fac;
            weight=ones(size(cube));
            logFlag=true; %cube=log10(cube);
            bartag='$\log n_{gas}\,[\mathrm{cm^{-3}}]$';
            printTypeTag='numberDensity';
        case {'dmdensity','rhodm','dm'}
            cube=smooth3(RHODM(boxx),'gaussian',7,3);
            
            weight=ones(size(cube));
            logFlag=true; %cube=log10(cube);
            bartag='$\log \rho_{DM}\,[\mathrm{M_\odot/Mpc^3}]$';
            printTypeTag='dmdensity';
        case {'totaldensity','rhotot','tot'}
            cube=RHOTOT(boxx);
            weight=ones(size(cube));
            logFlag=true; %cube=log10(cube);
            bartag='$\log \rho_{tot}\,[\mathrm{M_\odot/Mpc^3}]$';
            printTypeTag='totdensity';
        case {'entropy','ent','k','s'}
            %SVIR=TVIR.*RVIR.^2./(3.*MVIR./(4.*pi)).^(2./3); %normalization for entropy
            cube=S(boxx).*Units.factors.f_ent; % in units of KeV cm^2
            weight=rog;
            logFlag=true; %cube=log10(cube);
            bartag='$\log S\,[\mathrm{KeV\,cm^2}]$';
            printTypeTag='entropy';
        case {'temperature','temp','t'}
            cube=T(boxx);
            weight=rog;
            logFlag=true; %cube=log10(cube);
            bartag='$\log T\,[\mathrm{K}]$';
            printTypeTag='temperature';
        case {'mach'}
            
            if ~all(vcubeFlag)
                [hubX, hubY, hubZ] = hubble_flow(boxx,[0,0,0]);
                vx = Vx(boxx)+hubX-vcm(1);
                vy = Vy(boxx)+hubY-vcm(2);
                vz = Vz(boxx)+hubZ-vcm(3);
                vcubeFlag(:)=true;
            end
            vc=sqrt(vx.^2+vy.^2+vz.^2);
            
            
            rcube=mk_radius_cube(boxx); % r^2 cube in Mpc
            tc=read_T_Profile(rcube);
            %moo=0.5926;
            gm=5/3;
            cso=(gm*Units.kb/Units.mm.*tc).^0.5;  % ./1e5;  % so
            cso=cso./Units.km;
            cube=vc./cso;
            
            %logFlag=true; %cube=log10(cube);
            bartag='$\mathcal{M}$';
            %bartag='$\log Mach$';
            printTypeTag='mach';
            
        case {'pressure','press','p'}
            cube=T(boxx);
            cube=cube.*rog./(TVIR*MVIR/(4*pi*RVIR^3/3));
            logFlag=true; %cube=log10(cube);
            bartag='$\log P\,[\mathrm{P_{vir}}]$';
            printTypeTag='pressure';
            
        case('rps')
            fac=(Units.Ms.*1e-3)/(Units.Mpc*1e-2)^3*(Units.km*1e-2)^2*1e3; % units of mPa 
            cube=alfa.*rog.*(get_vvir.*Units.km).^2.*fac;
            
            
            %rpscube./(TVIR*MVIR/(4*pi*RVIR^3/3));
            logFlag=true; %cube=log10(cube);
            bartag='$\log P_{\mathrm{ram}}\,[\mathrm{mPa}]$';
            weight=ones(size(cube));
            printTypeTag='rps';
            
            
        case {'xray','lx'}
%             tb=T(boxx);
%             tb(tb<1e5)=0;
%             cube=rog.^2.*sqrt(tb);
            %cube=cube.*rog;%./(TVIR*MVIR/(4*pi*RVIR^3/3));
            vol=(boxx/NCELL/hub/(1+zred))^3;
            cube=cooling_cube(boxx,'cgs/Mpc').*vol;
            logFlag=true; %cube=log10(cube);
            bartag='$\log L_X\,[\mathrm{10^{42}erg/sec}]$';
            printTypeTag='xray';
            weight=ones(size(cube));
            slTypeDef='sum';
        case {'iametals','zia'}
            cube=ZIa(boxx);
            logFlag=true;
            bartag='$\log(ZIa)\,[\mathrm{Z_\odot}]$';
            printTypeTag='metalicityZIa';
        case {'iimetals','zii'}
            cube=ZII(boxx);
            bartag='$\log(ZII)\,[\mathrm{Z_\odot}]$';
            logFlag=true;
            printTypeTag='metalicityZII';
        case {'metalicity','metals','ztot'}
            cube=ZIa(boxx)+ZII(boxx);
            bartag='$\log(Z)\,[\mathrm{Z_\odot}]$';
            logFlag=true;
            printTypeTag='metalicity';
        case {'potential','potent'}
            [pdm, pgs , ~, ~]=poteng(boxx);
            cube=(pgs+pdm);
            %edot=ro.*pdot.*vol;
            clear pdm pgs ro
            bartag='$\phi\,[\mathrm{(km/sec)^2}]$';
            printTypeTag='potential_energy';
        case {'kinetic'}
            bartag='$E_{kin}\,[\mathrm{(km/sec)^2}]$';
            printTypeTag='kinetic';
        case {'poteng'}
            vol=(boxx/NCELL/hub/(1+zred))^3;
            [pdm, pgs , ~, ro]=poteng(boxx);
            cube=(pgs+pdm).*ro.*vol;
            %edot=ro.*pdot.*vol;
            cube=-1*log10(abs(cube));
            clear pdm pgs pdot ro vol
            bartag='$E_{pot}\,[\mathrm{M_{\odot}(km/sec)^2}]$';
            printTypeTag='potential';
        case {'total_energy'}
        case {'cooling'}
        case {'velocity','vel'}
            if any(~vcubeFlag)
                [hubX, hubY, hubZ] = hubble_flow(boxx,[0,0,0]);
                vx = Vx(boxx)+hubX-vcm(1);
                vy = Vy(boxx)+hubY-vcm(2);
                vz = Vz(boxx)+hubZ-vcm(3);
                vcubeFlag(:)=true;
            end
            cube=sqrt(vx.^2+vy.^2+vz.^2);
            weight=rog;
            %logFlag=true; %cube=log10(cube);
            bartag='$v\,[\mathrm{km/sec}]$';
            printTypeTag='velocity';
    case {'losvelocity','losvel'}
            if any(~vcubeFlag)
                [hubX, hubY, hubZ] = hubble_flow(boxx,[0,0,0]);
                vx = Vx(boxx)+hubX-vcm(1);
                vy = Vy(boxx)+hubY-vcm(2);
                vz = Vz(boxx)+hubZ-vcm(3);
                vcubeFlag(:)=true;
            end
            
            if(sum(plotproj))==1
                if plotproj(1)==1
                    cube=vx;
                elseif plotproj(2)==1
                    cube=vy;
                elseif plotproj(3)==1
                    cube=vz;
                else
                    error('mkmap: something wrong with projecitons')
                end
            else
                error('mkmap: only one projection at a time for los velocity')
            end
                       
            weight=rog;
            %logFlag=true; %cube=log10(cube);
            bartag='$v\,[\mathrm{km/sec}]$';
            printTypeTag='velocity';
       
        case {'vr','radialvel'}
            if any(~vcubeFlag)
                [hubX, hubY, hubZ] = hubble_flow(boxx,[0,0,0]);
                vx = Vx(boxx)+hubX-vcm(1);
                vy = Vy(boxx)+hubY-vcm(2);
                vz = Vz(boxx)+hubZ-vcm(3);
                vcubeFlag(:)=true;
            end
            
            [meshY, meshX, meshZ] = meshgrid(1:NCELL, 1:NCELL, 1:NCELL);
            %convert to center origin coordinates
            meshX = meshX - (NCELL+1)/2 -cm(1);
            meshY = meshY - (NCELL+1)/2 -cm(2);
            meshZ = meshZ - (NCELL+1)/2 -cm(3);
            % Fix Units (to be in Mpc)
            meshX = meshX * ((boxx/hub)/NCELL);
            meshY = meshY * ((boxx/hub)/NCELL);
            meshZ = meshZ * ((boxx/hub)/NCELL);
            
            norm = sqrt(meshX.^2 + meshY.^2 + meshZ.^2);
            meshX = meshX./norm;
            meshY = meshY./norm;
            meshZ = meshZ./norm;
            
            cube=vx.*meshX + vy.*meshY + vz.*meshZ;
            
            clear meshX meshY meshZ norm
            
            weight=rog;
            %logFlag=true; %cube=log10(cube);
            bartag='$v_r\,[\mathrm{km/sec}]$';
            printTypeTag='radvelocity';
            
            
            
            
            
        otherwise
            disp('unknown type');
            return
    end
end

cube=cube./normfactor;

if ~isempty(contrType)
    switch lower(contrType)
        case{'none','off'}
            contrFlag=false;
        case {'density','gasdensity','rho','gas','rhogas'}
            contrCube=rog;
            contrWeight=ones(size(contrCube));
            contrLogFlag=true;
        case {'dmdensity','rhodm','dm'}
            contrCube=smooth3(RHODM(boxx),'gaussian',9,5);
            
            contrWeight=ones(size(contrCube));
            contrLogFlag=true;
        case {'totaldensity','rhotot','tot'}
            contrCube=RHOTOT(boxx);
            contrWeight=ones(size(contrCube));
            contrLogFlag=true;
        case {'entropy','ent','k','s'}
            %SVIR=TVIR.*RVIR.^2./(3.*MVIR./(4.*pi)).^(2./3); %normalization for entropy
            contrCube=S(boxx).*f_ent; % in units of KeV cm^2
            contrWeight=rog;
            contrLogFlag=true;
        case {'temperature','temp','t'}
            contrCube=T(boxx);
            contrWeight=rog;
            contrLogFlag=true;
            
        case {'pressure','press','p'}
            contrCube=T(boxx).*rog./(TVIR*MVIR/(4*pi*RVIR^3/3));
            contrLogFlag=true;
            contrWeight=ones(size(contrCube));
            
        case {'xray'}
      
            %tb=T(boxx);
            %tb(tb<1e5)=0;
            vol=(boxx/NCELL/hub/(1+zred))^3;
            contrCube=cooling_cube(boxx,'cgs/Mpc').*vol; %rog.^2.*sqrt(tb);
            contrWeight=ones(size(cube));
            contrLogFlag=true; %cube=log10(cube);
            slTypeCont='sum';
        case {'metals','zia'}
            contrCube=ZIa(boxx);
            
        case {'iametals','zii'}
            contrCube=ZII(boxx);
            
        case {'metalicity','iimetals','ztot'}
            contrCube=ZIa(boxx)+ZII(boxx);
            
        otherwise
            error('MKMAP - Illegal contour tpye: %s',contrType);
    end
end


%% prepare map parameters

%define slice index
if thick>0
    thk=ceil(0.5.*thick./boxx.*NCELL);   %% thick is in comoving Mpc/h
else
    thk=0.5*abs(thick);   %% defualt value of 6 cells for slice
end

if contThick>0
    cthk=ceil(0.5.*contThick./boxx.*NCELL);   %% thick is in comoving Mpc/h
else
    cthk=thk;
end

slind=(floor(0.5.*NCELL)-thk+1):1:(ceil(0.5.*NCELL)+thk);
slindCont=(floor(0.5.*NCELL)-cthk+1):1:(ceil(0.5.*NCELL)+cthk);

% % if one wants to shift the center of the image
% if cshift~=0
%     shift_cen=0.5.*NCELL-(ceil(0.5.*NCELL*(cshift/boxx*2+1)));
% else
%     shift_cen=0;
% end
shift_cen=0;

slind=slind-shift_cen;

% Zoom-in picture zoomBox is in Mpc/h
%zoomBox=0;
% if zoomFlag
%     zmind=ceil(zoomBox./boxx.*NCELL);
% else
%     zmind=0.5.*NCELL;
% end
zmind=0.5*NCELL;


sidind=(0.5*NCELL+1)-zmind:0.5*NCELL+zmind;
actbox=zmind*(boxx/NCELL);
side=[-actbox actbox];
if(~comoveFlag)
    side=side./(1+zred);
    actbox=actbox./(1+zred);
end
zoomBoxlen=length(sidind);

if ~exist('slType','var')
    slType=slTypeDef;
end

 slice=mk_slice(cube,weight,slind,slType);
if logFlag
    slice=log10(slice);
end

%% prepare contours
if ~exist('slTypeCont','var')
    slTypeCont=slTypeDef;
end
% prepare contour slice
if contrFlag
    contrSlice=mk_slice(contrCube,contrWeight,slindCont,slTypeCont);
    if contrLogFlag
        contrSlice=log10(contrSlice);
    end
end
cBsize=diff(side)./NCELL;
cside=side(1)+0.5*cBsize:cBsize:side(2)-0.5*cBsize;

%% velocity field parameters (default filte=4)
if velFlag2
    [v_x, v_y]=mk_vfield(vx,vy,vz,rog,slind);
    clear vx vy vz;
    diluted_len = length(1:dilute:zoomBoxlen);
    diluted_jump = 2.*actbox/(diluted_len-1);
    notdiluted_jump = 2.*actbox/(zoomBoxlen-1);
    [xxv, yyv] = meshgrid(-actbox:diluted_jump:actbox, -actbox:diluted_jump:actbox);
    [xxs, yys] = meshgrid(-actbox:notdiluted_jump:actbox,-actbox:notdiluted_jump:actbox);
end




%% draw map
for projection = 1:3
    if plotproj(projection)
        figure;
        if ~clim_Flag
            clims(1) = min(slice(:));
            clims(2) = max(slice(:));
            
        end
        
        % fix to avoid problems with setting the colorbar
        if logFlag
            clims(isinf(clims))=-30;
        end
        
        %% start plotting
        
        imagesc(side,side,squeeze(slice(:,:,projection)),clims);%axis equal tight;
        hold on
        %% add velocity field
        if velFlag2
            resampled_v_x = imresize(squeeze(v_x(:,:,projection)), [diluted_len diluted_len], 'bilinear');
            resampled_v_y = imresize(squeeze(v_y(:,:,projection)), [diluted_len diluted_len], 'bilinear');
            if strmFlag
                h = streamslice(xxs,yys,squeeze(v_x(sidind,sidind,projection)),squeeze(v_y(sidind,sidind,projection)),streamDense);
                set(h,'color',streamColor)
            end
            if velFlag
                quiver(xxv,yyv,resampled_v_x,resampled_v_y, 1.3,'color',vfieldColor);
            end
        end
        
        %% add contours
        
        if contrFlag
            scontrSlice=squeeze(contrSlice(:,:,projection));
            
            if isempty(contLevs)
                contLevs=linspace(min(min(scontrSlice)),max(max(scontrSlice)),contSpacing);
            elseif length(contLevs)==1
                contLevs=linspace(min(min(scontrSlice)),max(max(scontrSlice)),contLevs);
            end
            
            [~,ch]=contour(cside,cside,scontrSlice,contLevs,contColor,'linewidth',2);
            if contShowLabel
                set(ch,'ShowText','on','TextStep',get(ch,'LevelStep')*2)
            end
            %clabel(C,chh,'fontsize',18,'rotation',0)
        end
        
        % make sure marked circles fit in box
        if ~isempty(marks)
            if(comoveFlag)
                markss=marks.*(1+zred);
            else
                markss=marks;
            end
            marks_inbox=markss; %(marks<=(0.5.*boxx./hub));
            draw_circles(marks_inbox,markColor);
        end
        
        rv(1)=get_rvir(500);
        rv(2)=get_rvir();
        rv(3)=get_rvir(200);
        
        if(comoveFlag)
            rv=rv.*(1+zred);
        end
        draw_circle_boxx(gcf,rv(rvin),'white');      % draw rvir
        
%         %% draw circle 
%          if exist('circleStruct','var')
%             for indL=1:length(circleStruct) % go over all lines
%                 
%                 if ~isfield(circleStruct(indL),'color') % set color  if none is given
%                     circleStruct(indL).color='k';
%                 end
%                 if ~isfield(circleStruct(indL),'width') % set linewidth  if none is given
%                     circleStruct(indL).width=1.5;
%                 end
%                 if ~isfield(circleStruct(indL),'type') % set line type  if none is given
%                     circleStruct(indL).type='-';
%                 end
%                 circleVal=circleStruct(indL).value;
%                
%                 %draw line
%                 if(comoveFlag)
%                     circleVal=circleVal.*(1+zred);
%                 end
%                 draw_circle_boxx(gcf,circleVal,circleStruct(indL).color);      % draw rvir
%                 
%             end
%             
%          end
         
        
        
        
        
        
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
            lunit='[\mathrm{Mpc/h\ comoving\,}]';
        else
            lunit='[\mathrm{Mpc/h}]';
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
        
        %% sort out ticks
        if ~exist('tickStruct','var')
            tickjump=0.25.*boxx;
            
            if ~comoveFlag && strcmp(aexpn,'a06')
                tickjump=0.15.*boxx;
            end
            ticks=-2*tickjump:tickjump:2*tickjump;
            if ~exist('xticks','var')
                xticks=ticks;
            end
            if ~exist('yticks','var')
                yticks=ticks;
            end
            if boxx==1
                xtickLabels=sprintf('%3.2f|',xticks);
                ytickLabels=sprintf('%3.2f|',yticks);
                zl='|0.00|';
                
            else  %if any(boxx==[2 4])
                xtickLabels=sprintf('%3.1f|',xticks);
                ytickLabels=sprintf('%3.1f|',yticks);
                zl='|0.0|';
                %             else
                %                 xtickLabels=sprintf('%3.0f|',xticks);
                %                 ytickLabels=sprintf('%3.0f|',yticks);
                %                 zl='|0|';
                %
            end
            
            % fix zero in labels
            k=strfind(xtickLabels,zl);
            if k>0
                xtickLabels=sprintf('%s%s%s',xtickLabels(1:k-1),'|0|',xtickLabels(k+length(zl):end));
            end
            
            k=strfind(ytickLabels,zl);
            if k>0
                ytickLabels=sprintf('%s%s%s',ytickLabels(1:k-1),'|0|',ytickLabels(k+length(zl):end));
            end
            
        else
            xticks=tickStruct.xticks;
            yticks=tickStruct.yticks;
            %xtickLabels=tickStruct.xtickLabels;
            %ytickLabels=tickStruct.ytickLabels;
        end
        
        
        
        set(gca,'Ydir','normal','Fontsize',14,...
            'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1],...
            'XTick',xticks ,'YTick',yticks,'TickLength',[-0.015 -0.015],...
            'XLim',side,'YLim',side);    %'XTickLabel',xtickLabels,'YTickLabel',ytickLabels
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
        %else
        %    set(gca,'XLim',[-boxx./2 boxx./2],'YLim',[-boxx./2 boxx./2])
        %end
        
        
        %   if mvirboxFlag
        %       ss=mk_mvir_string(MVIR);
        %       text(strpos(1),strpos(2),ss,'Interpreter','latex','fontsize',14);
        %   end
        
        %% do bar and colormap stuff
        caxis(clims);
        
        if brewerFlag
            map = brewermap(256,brewMap);
        end
        colormap(map);  
        %colormap(avijet);  % FIX? 
        bar=colorbar;
        if ~isempty(bartagOveride) || ~exist('bartag','var')
            bartag=bartagOveride;
        end
        barTitle(bar,bartag);
        set(bar,'Fontsize',12)
        %set(get(bar,'Title'),'String',bartag,'Fontsize',12,'Interpreter','latex');
        set(gcf,'Colormap',map);
        %set(gcf,'Colormap',avijet);
        %title(sprintf('%s %s, Thickness=%s Mpc/h',CLUSTER,type,num2str(thick,3)),'Fontsize',12,'Interpreter','latex');
        
        %% set title
        if titleFlag
            switch titletype
                case 'full'
                    mvstr=mk_mvir_string(MVIR);
                    titlemine(sprintf('%s %s, %s $z=%3.2g$',CLUSTER,mvstr,titletag,zred));
                case 'custom'
                    titlemine(titletag);
                otherwise
                    error('MKMAP: Illegal value for title type (full/custom)')
            end
        end
        
        %% add grid
        if gridFlag
            grid;
        end
        
        %set(gcf,'Position',[ 560   355   727   593]);
        
        %set(gcf,'Colormap',avijet);%'PaperOrientation','landscape')
        
        %printing:
        %Don't want more than 3 printing arugments altogether
        %numvarargs=length(varargin);
        %if numvarargs>3
        %    error('mkmap: too many printing arguments');
        %end
        
        %set defualt printing Flag, printing tag, and output dir
        %defvals={'noprint', type ,'/home/eladzing/Ubuntu One/cluster/printout'};
        
        %assign the optional values
        %defvals(1:numvarargs)=varargin;
        
        % transfer to easy to use varaibles
        %[printFlag printtag printoutdir]=defvals{:};
        
        %% print
        if printFlag
            name=sprintf('%s_map%s_b%d_%s_%s_%s',CLUSTER,prjtag,boxx,printTypeTag,printtag,aexpn);
            %name=sprintf('%s/%s_map%s_b%d_%s.%s',printoutdir,CLUSTER,prjtag,boxx,printtag,'%s');
            
            switch(printformat)
                case 'pdf'
                    printout_fig(gcf,name,'dir',printoutdir,'v')
                case 'png'
                    printout_fig(gcf,name,'dir',printoutdir,'v','png')
                otherwise
                    error('mkmap - Illegal print format %s',printformat)
            end
            
            %             if printformat(1)
            %                 exportfig(gcf,sprintf(name,'png'),'format','png');
            %             end
            %             if printformat(2)
            %                 exportfig(gcf,sprintf(name,'eps'));
            %             end
            %saveas(gcf,sprintf('%s/%s_map%s_b%d_%s.png',printoutdir,CLUSTER,prjtag,boxx,ptag));
        end
        
    end
end
%
% 