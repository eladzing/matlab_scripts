function mkmap_rps(boxx,varargin)

% basic function for making a map of a cluster plus a velocity field.
% boxx -which boxx to map
% plotproj - which projection to plot, logical vector
% type - which value do we plot
%clims - limits of colorbar
%velflag - flag to  velocity field or not.
% vcm - center of mass velocity should be calculated in advance
% thick - the thickness of the slice in Mpc/h comoving
% marks - a vector of radii for drawing circles in units of Mpc
% dilute - vector field dilution. a value <0 gives the defualt dilute=4

%% preamble

% argument test
if nargin < 4
    error('mkmap: not enough arguments')
    %elseif nargin > 10
    %    error('mkmap: too many arguments')
end

units;

global CLUSTER;
global VCM
%global HALO_PATH;
global NCELL;
global zred
global hub;
global DEFUALT_PRINTOUT_DIR
global aexpn;
%load(sprintf('%s/virial%d_%s', HALO_PATH, 8,aexpn));
projectiontags={'YZ' 'XZ' 'XY'};
bartagOveride='';

RVIR=get_rvir();
MVIR=get_mvir();
VVIR=get_vvir();
TVIR=get_tvir();

% defualt values
%boxx=8;
plotproj=[1 1 1]; % plot all three projections
dilute=-1;
marks=[];
vcm=[];
printoutdir=DEFUALT_PRINTOUT_DIR;
printtag='badtag';
printformat=[1 1];
type='';
thick=-6   ;
normfactor=1;
titletag='';
labeltype='full';
titletype='full';
cont_type='none';
cont_vals=[];
rvin=logical([0 1 0]); % must be 3 component vector of 0 or 1 for plotting r500 rvir r200.
markColor='k';
% default is to plot rvir alone.

% rps values
stripmodel='mv_strip'; %stripping model shown - default is the stripping radius for a given sat mass.
rps_vel='vvir';
fgas=0.1;  % baryon fraction within halo
alfa=1.0; %rps efficiency parameter
mvref=1e12;
cont_param=[];

%pval=-6;

%flags
cubeflag=false;
typeflag=false;
wtflag=false;
%labelflag=false;
titleflag=true;
printflag=false;
clim_flag=false; % no colorbar limits given, use min max
velflag=logical([0 0]); % first index controls arrows, second controls streamlines
%mvirboxflag=false;
logflag=false;
comoveflag=true;
contour_flag=false;
rpscubeflag=false;
vcubeflag(1:3)=false;
rogflag=false;
verboseFlag=false;
brewerFlag=false;

load('MyColormaps','newJet');
map=newJet;
%% reading arguments
i=1;
while i<=length(varargin)
     if verboseFlag  %show next argument 
        disp(varargin{i})
    end
    switch lower(varargin{i})
        case {'verbose','v'}
            verboseFlag=true;
%         case {'boxx','box'}   %which box to plot
%             i=i+1;
%             boxx=varargin{i};
        case{'data','datacube'} %plot a given cube;
            i=i+1;
            cubeflag=true;
            cube=varargin{i};
        case{'rpscube'} % rps cube essential for calculating contours
            i=i+1;
            rpscubeflag=true;
            rpscube=varargin{i};
        case{'rocube','rog'} % rps cube essential for calculating contours
            i=i+1;
            rogflag=true;
            rog=varargin{i};
        case{'weightcube','wtcube','wtdatacube','wt'} % weight for cube
            i=i+1;
            wtflag=true;
            weight=varargin{i};
        case{'datatype','type'} %prepare data cube of this type
            i=i+1;
            typeflag=true;
            type=varargin{i};
        case{'vxcube'};
            i=i+1;
            vcubeflag(1)=true;
            vx=varargin{i};
        case{'vycube'};
            i=i+1;
            vcubeflag(2)=true;
            vy=varargin{i};
        case{'vzcube'};
            i=i+1;
            vcubeflag(3)=true;
            vz=varargin{i};
        case{'contourtype','contour'} %Which contorur to draw (if any)
            i=i+1;
            contour_flag=true;
            cont_type=varargin{i};
        case{'contour_values','contvals'} %Which contorur to draw (if any)
            i=i+1;
            cont_vals=varargin{i};
        case{'contour_paramter','cont_param'}
            i=i+1;
            cont_param=varargin{i};
        case{'stripping_model','strip_model'} %type of conture - mass or stripping radius
            i=i+1;
            stripmodel=varargin{i};
        case{'strip_parameter','strip_param'}
            i=i+1;
            stripmodel_param=varargin{i};
        case{'alfa','efficency','alpha'} %rps efficiency parameter
            i=i+1;
            alfa=varargin{i};
        case{'model','satmodel','sat'}
            i=i+1;
            model=varargin{i};
        case{'rps_velocity','rps_v'}
            i=i+1;
            rps_vel=varargin{i};
        case{'gasfrac','fgas'}
            i=i+1;
            fgas=varargin{i};
        case{'log','log10'}
            logflag=true;
        case{'width','thick'} %thicknes of slice in Mpc/h comoving (like box)
            i=i+1;
            thick=varargin{i};
        case {'proj','projection'}
            i=i+1;
            plotproj=varargin{i};
            case {'yz','zy'}
            plotproj(1)=1;
        case {'zx','xz'}
            plotproj(2)=1;
        case {'xy','yx'}
            plotproj(3)=1;
        case 'all'
        case {'clims','limits'}
            i=i+1;
            clim_flag=true;
            clims=varargin{i};
            if diff(clims)==0
                clim_flag=false;
            end
        case {'velflag','velocity','vfield'}
            velflag=true;
            %             i=i+1;
%             answer=varargin{i};
%             if ischar(answer)
%                 velflag(:)=(strcmp(answer,'yes') || strcmp(answer,'on'));
%             else
%                 velflag(:)=logical(answer);
%             end
            
        case {'dilute','arrowspace'}
            i=i+1;
            dilute=varargin{i};
        case {'vcm'}
            i=i+1;
            vcm=varargin{i};
        case {'normalize','norm'}
            i=i+1;
            normfactor=varargin{i};
        case {'marks','circles'}
            i=i+1;
            marks=varargin{i};
        case {'rv','rv_circles'}   % must be 3 component vector of 0 or 1 for plotting r500 rvir r200.
            i=i+1;
            rvin=logical(varargin{i});
               
        case {'bartag','bartitle'}
            i=i+1;
            bartagOveride=varargin{i};
        case {'labels'}
            i=i+1;
            labeltype=varargin{i};
            %labelflag=(strcmp(l,'yes') || strcmp(answer,'on') || strcmp(answer));
        case {'title'}
            i=i+1;
            titletag=varargin{i};
            titleflag=~(strcmp(titletag,'no') || strcmp(titletag,'off'));
        case {'titletype'}
            i=i+1;
            titletype=varargin{i};
        case {'bartag','bartitle'}
            i=i+1;
            bartag=varargin{i};
            %          case {'mvirbox'} % place a text of Mvir on the map
            %             mvirboxflag=true;
            %             bartag=varargin{i+1};
            %          case {'mvirbox'} % place a text of Mvir on the map
            %             mvirboxflag=true;
            %             strpos=varargin{i+1};  % position of tex            strpos=varargin{i+1};  % position of text in map coordinates
        case {'comoving','comove','comoving_coord'}
            comoveflag=true;
        case {'proper','proper_coord'}
            comoveflag=false;
        case {'brewer','brewermap'}
            i=i+1;
            brewMap=varargin{i};
            brewerFlag=true;
        case {'print'}
            printflag=true;
        case {'outputdir','printout','printoutdir'}
            i=i+1;
            printoutdir=varargin{i};
        case {'format','printformat'}
            i=i+1;
            answer=varargin{i};
            switch answer
                case {'png','PNG'}
                    printformat=[1 0];
                case {'EPS','eps'}
                    printformat=[0 1];
                case {'both'}
                    printformat=[1 1];
                otherwise
                    error('mkmap printformat: illegal argument %s',varargin{i})
            end
            
        case {'printtag','printag'}
            i=i+1;
            printtag=varargin{i};
       otherwise
            error('mkmap_rps:illegal argument %s',varargin{i})
    end
    i=i+1;
   
   
end

%% argument tests
if ~any(boxx==[1 2 4 8])
    error('mkmap_rps:illegal box')
end

if ~(cubeflag || typeflag)
    error('mkmap: must enter data or datatype')
end
if (cubeflag && typeflag)
    error('mkmap: too many data  arguments')
end

if (length(vcm)~=3)
    vcm = VCM;
end

if ~wtflag
    weight=ones(NCELL,NCELL,NCELL);
end

if (strcmp(type,'mv_strip') || strcmp(type,'mv_strip_inv'))
    if ~strcmp(stripmodel,'mv_strip')
        error('mkmap_rps: incompatible stripping model for data type')
    end
elseif strcmp(type,'rstrip')
    if ~strcmp(stripmodel,'rstrip')
        error('mkmap_rps: incompatible stripping model for data type')
    end
end

if (isempty(cont_param) && contour_flag)
    if (strcmp(stripmodel,cont_type) || (strcmp(stripmodel,'mv_strip') && strcmp(cont_type,'mv_strip_inv')))
        cont_param=stripmodel_param;
    else
        error('nknmap_rps: no parameter for contour is given')
    end
end
%% begin preparing data
%bcson=sqrt(1.52e8.*T(boxx))./1e5; %speed of sound in km/sec

cm=[0,0,0];

if(any(velflag) || strcmpi(type,'flux'))
    [hubX hubY hubZ] = hubble_flow(boxx,[0,0,0]);
    if ~all(vcubeflag)
        vx = Vx(boxx);
        vy = Vy(boxx);
        vz = Vz(boxx);
    end
    
    vx = vx+hubX-vcm(1);
    vy = vy+hubY-vcm(2);
    vz = vz+hubZ-vcm(3);
    if ~rogflag
        rog=RHOG(boxx);
        rogflag=true;
    end
end

% prepare rps cube
if ~rpscubeflag
    if ~rogflag
        rog=RHOG(boxx);
        rogflag=true;
    end
    switch lower(rps_vel)
        case('vvir')  % assume sat moves at vvir
            rpscube=rog.*VVIR.^2;
        case('circ')  % assume sat moves at circular velocity G.*M/r
            %evalute v^2 as GM(<r)/r
            [meshY, meshX, meshZ] = meshgrid(1:NCELL, 1:NCELL, 1:NCELL);
            %convert to center origin coordinates
            meshX = meshX - (NCELL+1)/2 -cm(1);
            meshY = meshY - (NCELL+1)/2 -cm(2);
            meshZ = meshZ - (NCELL+1)/2 -cm(3);
            % Fix Units (to be in Mpc)
            meshX = meshX * ((boxx/hub)/NCELL);
            meshY = meshY * ((boxx/hub)/NCELL);
            meshZ = meshZ * ((boxx/hub)/NCELL);
            
            rcube=sqrt(meshX.^2+meshY.^2+meshZ.^2); % r cube in Mpc
            mtot=read_MTOT_Profile(rcube);
            fac=G;
            rpscube=rog.*mtot./rcube.*fac;
            clear rcube mtot
            
        otherwise
            error('mkmap_rps: illeagal velocity type for rps');
    end
end

%set up stripping model

switch stripmodel
    %calculate the stripping radius for a given sat mv
    case {'mv_strip','mv_strip_inv'}
        % get virial data for sat
        [rv mv , ~, ~]=calculate_virials('mvir',stripmodel_param);
        
        pval=(fgas/alfa)*(G*Ms/Mpc/km^2)/(4*pi)*(mv^2/rv^4);
        
    case 'rstrip'
        % calculate a given stripping radius for different sat masses
        % find vir constant
        
        mm=mvref;
        [rv ~ , ~, ~]=calculate_virials('mvir',mm);
        vir_const=mm^(1/3)/rv;
        
        pval=(fgas/alfa)*(G*Ms/Mpc/km^2)/(4*pi)*vir_const^4/(stripmodel_param^2);
        
        %     case 'disk'
        %         %default values
        %         %sigma=4e7; fd=0.2; beta=1;
        %         md=1e11; rd= 20; fd=0.25;
        %
        %         siz=length(stripmodel_param);
        %         switch siz
        %             case 1
        %                 sigma=stripmodel_param(1); % stellar sigma = Md/(2*pi*rd^2)
        %             case 2
        %                 sigma=stripmodel_param(1);
        %                 fd=stripmodel_param(2);  % gas fraction in disk
        %
        %             case 3
        %                 sigma=stripmodel_param(1);
        %                 fd=stripmodel_param(2);
        %
        %                 beta=stripmodel_param(3);  % star to gas scale radius ratio
        %             otherwise
        %                 error('mkmap_rps: disk stripping must recieve vector of no more than 3 values');
        %         end
        %         pval=(2.0*G/pi*(fd/alfa)*(md*Ms)^2/(rd*kpc)^4)*(Mpc^3/Ms/km^2);
        %
    case {'expdisk' 'expdisk_sigma'}
        
        switch stripmodel
            case 'expdisk'
                %default values
                md=1e10; rd= 4; fd=0.1; beta=l;
                
                switch length(stripmodel_param);
                    case 1
                        error('mkmap_rps: must enter at least disk mass & scale radius');
                    case 2
                        md=stripmodel_param(1);
                        rd=stripmodel_param(2);
                    case 3
                        md=stripmodel_param(1);
                        rd=stripmodel_param(2);
                        fd=stripmodel_param(3);
                        
                    case 4
                        md=stripmodel_param(1);
                        rd=stripmodel_param(2);
                        fd=stripmodel_param(3);
                        beta=stripmodel_param(4);
                    otherwise
                        error('mkmap_rps: disk stripping must recieve vector of no more than 3 values');
                end
                
                sigma=(md.*Ms)./(2*pi*(rd.*kpc)^2);
                
            case 'expdisk_sigma'
                %default values
                sigma=4e7; fd=0.1; beta=1;
                
                switch length(stripmodel_param)
                    case 1
                        sigma=stripmodel_param(1); %sigma_effective in Msun/kpc^2
                    case 2
                        sigma=stripmodel_param(1); %sigma_effective in Msun/kpc^2
                        fd=stripmodel_param(2); %rd_gas/rd_star
                    case 3
                        sigma=stripmodel_param(1); %sigma_effective
                        fd=stripmodel_param(2); %rd_gas/rd_starin Msun/kpc^2
                        beta=stripmodel_param(3); %gas fraction
                    otherwise
                        error('mkmap_rps: disk stripping must recieve vector of mass and radius');
                end
        end
        
        sigma=sigma.*(Ms/kpc^2);
        
        % create force vs strip calculation
        r=0.01:0.01:60;
        mf1=exp_disk_mass(r,beta);
        fg1=disk_force_reduced(r,'beta',beta,'fg',0.1);
        [fgmax,ind1]=max(fg1);
        fg=fg1(ind1:end);
        mf=mf1(ind1:end);
        fg_mf01=interp1(mf(mf<0.5),fg(mf<0.5),0.1); % value of fg for mf=0.1
        
        clear r ind1 mf1 fg1
        
        pval0=pi*G*sigma^2*fd*(Mpc^3/Ms/km^2); % prefactors of the disk force per unit area
        pvalS=pi*G*fd.*fg_mf01;
        %pval=G*pi*fd/alfa*(re_fac^2*sig*Ms/kpc^2)^2/beta^2*(Mpc^3/Ms/km^2);
        %pval=(2.0*G/pi*(fd/alfa)*(md*Ms)^2/(rd*kpc)^4)*(Mpc^3/Ms/km^2);
        
        
    otherwise
        error('mkmap_rps: illegal value for stripping model');
end


if typeflag
    if ~rogflag
        rog=RHOG(boxx);
        rogflag=true;
    end
    switch lower(type)
        case('mv_strip') % show stripping radius for a given sat mass
            cube=1./sqrt(rpscube./pval);
            logflag=false;
            bartag='$\frac{R_{strip}}{R_{sat}}$';
            weight=ones(size(cube));
        case('mv_strip_inv') % show stripping radius for a given sat mass
            cube=1./sqrt(rpscube./pval);
            cube(cube>1)=1;
            cube=100.*(1-cube);
            lc=-5:10:100;
            logflag=false;
            bartag='Strip $\%$';
            weight=ones(size(cube));
            bticks=0:10:100;
             printTypeTag='rpsM';
        case('rstrip') % given stripping radius for different sat masses
            cube=(rpscube./pval).^(3/2);
            logflag=true;
            lc=6:0.5:12;
            bticks=8:0.5:12;
            bartag='$M_{sat}\,[\mathrm{M_\odot}]$';
            weight=ones(size(cube));
            printTypeTag='rpsR';
            %         case('rps_disk')
            %             cube=rpscube./pval;
            %             logflag=true;
            %             bartag='$\log(\zeta_{disk})$';
            %             weight=ones(size(cube));
            
        case('expdisk')
            pv=alfa.*rpscube./pval0;
            cube=interp1(fg,mf,pv);
            cube(pv>fgmax)=0.0;
            cube=(1-cube).*100;
            logflag=false;
            bartag='$\%M_{strip}$';
            weight=ones(size(cube));
            clear pv fg mf
            
        case('expdisk_sig')
            cube=sqrt(alfa.*rpscube./pvalS);
            logflag=true;
            bartag='$\log(\Sigma_s)$';
            weight=ones(size(cube));
            clear pv
            
        case('expdisk_sig_eff')
            cube=sqrt(alfa.*rpscube./pvalS).*sigma_factor('half');
            logflag=true;
            bartag='$\log(\Sigma_{\mathrm{eff})$';
            weight=ones(size(cube));
            clear pv
            
        case('expdisk_sig_5090')
            cube=sqrt(alfa.*rpscube./pvalS).*sigma_factor('5090');
            logflag=true;
            bartag='$\log(\Sigma_{5090})$';
            weight=ones(size(cube));
            clear pv
            
            
        case('rps')
            cube=rpscube./(TVIR*MVIR/(4*pi*RVIR^3/3));
            logflag=true; %cube=log10(cube);
            bartag='$\log P_{ram}\,[\mathrm{P_{vir}}]$';
            weight=ones(size(cube));
        case('rps_over_pressure')
            cube=rpscube./(T(boxx).*rog.*f_pr);
            logflag=true; %cube=log10(cube);
            bartag='$\log \frac{P_{ram}}{P_{ICM}} $';
            weight=ones(size(cube));
        case {'flux','massflux'}
            fluxnorm=(0.056.*(MVIR./1e13).^0.15)./(4.*pi); % normalization for flux
            
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
            [Mg200 , ~, ~]=read_Mass_Profiles(RVIR);
            cube=rog.*vr.*rcube2./Mg200.*(1e5./3.0856e24.*3.15e7.*1e9)./fluxnorm ;
            weight=ones(size(cube));
            
            bartag='$\frac{\dot M/M_{gas}}{d\Omega} [\frac{\dot M}{M_{gas}}\left|_{\mathrm{virial}}]$';
            
            clear vr rcube2 meshX meshY meshZ
            
        case {'density','gasdensity','rho','gas','rhogas'}
            cube=rog;
            
            weight=ones(size(cube));
            logflag=true; %cube=log10(cube);
            bartag='$\log \rho_{gas}\,[\mathrm{M_\odot/Mpc^3}]$';
        case {'dmdensity','rhodm','dm'}
            cube=RHODM(boxx);
            weight=ones(size(cube));
            logflag=true; %cube=log10(cube);
            bartag='$\log \rho_{DM}\,[\mathrm{M_\odot/Mpc^3}]$';
        case {'totaldensity','rhotot','tot'}
            cube=RHOTOT(boxx);
            weight=ones(size(cube));
            logflag=true; %cube=log10(cube);
            bartag='$\log \rho_{tot}\,[\mathrm{M_\odot/Mpc^3}]$';
        case {'entropy','ent','k','s'}
            SVIR=TVIR.*RVIR.^2./(3.*MVIR./(4.*pi)).^(2./3); %normalization for entropy
            cube=S(boxx)./SVIR;
            weight=rog;
            logflag=true; %cube=log10(cube);
            bartag='$\log K [\mathrm{virial}]$';
        case {'temperature','temp','t'}
            cube=T(boxx);
            weight=rog;
            logflag=true; %cube=log10(cube);
            bartag='$\log T\,[\mathrm{K}]$';
        case {'pressure','press','p'}
            cube=T(boxx);
            cube=cube.*rog./(TVIR*MVIR/(4*pi*RVIR^3/3));
            logflag=true; %cube=log10(cube);
            bartag='$\log P\,[\mathrm{P_{vir}}]$';
        case {'metals','zia'}
            cube=ZIa(boxx);
            bartag='ZIa';
        case {'iametals','zii'}
            cube=ZII(boxx);
            bartag='ZII';
        case {'metalicity','iimetals','ztot'}
            cube=ZIa(boxx)+ZII(boxx);
            bartag='Z';
        case {'potential','potent'}
            [pdm pgs , ~, ~]=poteng(boxx);
            cube=(pgs+pdm);
            %edot=ro.*pdot.*vol;
            clear pdm pgs ro
            bartag='$\phi\,[\mathrm{(km/sec)^2}]$';
        case {'kinetic'}
            bartag='$E_{kin}\,[\mathrm{(km/sec)^2}]$';
            
        case {'poteng'}
            vol=(boxx./NCELL./hub).^3;
            
            [pdm pgs , ~, ro]=poteng(boxx);
            cube=(pgs+pdm).*ro.*vol;
            %edot=ro.*pdot.*vol;
            cube=-1*log10(abs(cube));
            clear pdm pgs pdot ro vol
            bartag='$E_{pot}\,[\mathrm{M_{\odot}(km/sec)^2}]$';
        case {'total_energy'}
        case {'cooling'}
        otherwise
            error('mkmap_rps: unknown type');
            
    end
end

cube=cube./normfactor;

%% prepare map parameters

%define slice index
if thick>0
    thk=ceil(0.5.*thick./boxx.*NCELL);   %% thick is in comoving Mpc
else
    thk=0.5*abs(thick);   %% defualt value of 6 cells for slice
end
slind=(floor(0.5.*NCELL)-thk+1):1:(ceil(0.5.*NCELL)+thk);

% % if one wants to shift the center of the image
% if cshift~=0
%     shift_cen=0.5.*NCELL-(ceil(0.5.*NCELL*(cshift/boxx*2+1)));
% else
%     shift_cen=0;
% end
shift_cen=0;

slind=slind-shift_cen;
% Zoom-in picture zmbox is in Mpc/h
zmbox=0;
if zmbox~=0
    zmind=ceil(zmbox./boxx.*NCELL);
else
    zmind=0.5.*NCELL;
end

sidind=(0.5*NCELL+1)-zmind:0.5*NCELL+zmind;
actbox=zmind*(boxx/NCELL);
side=[-actbox actbox];
if(~comoveflag)
    side=side./(1+zred);
    actbox=actbox./(1+zred);
end
cside=side(1):diff(side)/(NCELL-1):side(2);
zoomboxlen=length(sidind);

slice=mk_slice(cube,weight,slind);
rps_slice=mk_slice(rpscube,weight,slind);

if logflag
    slice=log10(slice);
end

% velocity field parameters
if dilute<=0
    dilute = 4;
end
if any(velflag)
    [v_x v_y]=mk_vfield(vx,vy,vz,rog,slind);
    clear vx vy vz;
    diluted_len = length(1:dilute:zoomboxlen);
    diluted_jump = 2.*actbox/(diluted_len-1);
    notdiluted_jump = 2.*actbox/(zoomboxlen-1);
    [xxv yyv] = meshgrid(-actbox:diluted_jump:actbox, -actbox:diluted_jump:actbox);
    [xxs yys] = meshgrid(-actbox:notdiluted_jump:actbox,-actbox:notdiluted_jump:actbox);
end


tickjump = actbox/4;
load('MyColormaps','avijet');



% draw map
for projection = 1:3
    if plotproj(projection)
        figure;
        if ~clim_flag
            clims(1) = min(slice(:));
            clims(2) = max(slice(:));
            
        end
        
        % fix to avoid problems with setting the colorbar
        if logflag
            clims(isinf(clims))=-30;
        end
        
        hold on
        imagesc(side,side,squeeze(slice(:,:,projection)),clims);%axis equal tight;
        %ss=linspace(side(1),side(2),NCELL);
        %contourf(ss,ss,squeeze(slice(:,:,projection)),lc)
        if velflag(1)
            resampled_v_x = imresize(squeeze(v_x(:,:,projection)), [diluted_len diluted_len], 'bilinear');
            resampled_v_y = imresize(squeeze(v_y(:,:,projection)), [diluted_len diluted_len], 'bilinear');
            quiver(xxv,yyv,resampled_v_x,resampled_v_y, 2, 'k');
        end
%         if velflag(1)
%             h = streamslice(xxs,yys,squeeze(v_x(sidind,sidind,projection)),squeeze(v_y(sidind,sidind,projection)),1.5);
%             set(h,'color','k')
%         end
        
        %% draw contours
        if contour_flag
            cont_slice=squeeze(rps_slice(:,:,projection));
            switch cont_type
                case 'mv_strip'
                    [trv tmv , ~, ~]=calculate_virials('mvir',cont_param);
                    pv=(fgas/alfa)*(G*Ms/Mpc/km^2)/(4*pi)*(tmv^2/trv^4);
                    
                    rstr=1./(sqrt(cont_slice./pv));
                    conts_base=[0.01 0.05 0.1 0.2 0.5 0.75 1.0];
                    clear trv tmv pv
                case 'mv_strip_inv'
                    [trv tmv , ~, ~]=calculate_virials('mvir',cont_param);
                    pv=(fgas/alfa)*(G*Ms/Mpc/km^2)/(4*pi)*(tmv^2/trv^4);
                    
                    rstr=1./(sqrt(cont_slice./pv));
                    rstr(rstr>1)=1;
                    rstr=100.*(1-rstr);
                    conts_base=[1 5 10 25 50 75 90 95];
                    clear trv tmv pv
                case 'rstrip'
                    tmv=mvref;
                    [trv ~ , ~, ~]=calculate_virials('mvir',tmv);
                    vir_const=tmv^(1/3)/trv;
                    
                    pv=(fgas/alfa)*(G*Ms/Mpc/km^2)/(8*pi)*vir_const^4/(cont_param^2);
                    
                    rstr=log10((cont_slice./pv).^(3/2));
                    conts_base=[9 10 11 12 13];
                    clear trv tmv pv vir_const
                    
                case {'alpha','alfa'}
                    conts_base=[0.05 0.1 0.2 0.5 1.0 2];
                    switch stripmodel
                        case 'mv_strip'
                            rstr=(alfa*pval)./(cont_slice.*cont_param^2);
                        case 'rstrip'
                            rstr=((alfa*pval)*cont_param^(2/3))./(cont_slice);
                    end
                    
                otherwise
                    error('mkmap_rps: illegal value for contour type');
            end
            
            if isempty(cont_vals)
                conts=conts_base;
            else
                conts=cont_vals;
            end
            
            C=contour(cside,cside,rstr,conts,'-k','linewidth',1.5);
            clabel(C);
        end
        
        
        %% draw marks and circles
        if ~isempty(marks)
            if(comoveflag)
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
        
        if(comoveflag)
            rv=rv.*(1+zred);
        end
        draw_circle_boxx(gcf,rv(rvin),'white');      % draw rvir
        
       
        hold off
        
        %% units labels and titles
        prjtag=projectiontags{projection};
        if(comoveflag && zred>0)
            lunit='[\mathrm{Mpc/h\,comoving}]';
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
        
        tickjump=0.25.*boxx;
        
        if ~comoveflag && strcmp(aexpn,'a06')
            tickjump=0.15.*boxx;
        end
        ticks=-2*tickjump:tickjump:2*tickjump;
        xticks=ticks;
        yticks=ticks;
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
        
        
        set(gca,'Ydir','normal','Fontsize',14,...
            'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1],...
            'XTick',xticks ,'YTick',yticks,'TickLength',[-0.015 -0.015],...
            'XLim',side,'YLim',side);    %'XTickLabel',xtickLabels,'YTickLabel',ytickLabels
        box on;
        
        %   if mvirboxflag
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
        set(bar,'Fontsize',12,'ticks',bticks)
        %set(get(bar,'Title'),'String',bartag,'Fontsize',12,'Interpreter','latex');
        set(gcf,'Colormap',map);
        
        %% title
        if titleflag
            switch titletype
                case 'full'
                    mvstr=mk_mvir_string(MVIR);
                    titlemine(sprintf('%s %s, %s, $z=%3.2g$',CLUSTER,mvstr,titletag,zred));
                case 'half'
                    titlemine(sprintf('%s, %s, $z=%3.2g$',CLUSTER,titletag,zred));
                case 'custom'
                    titlemine(titletag);
                otherwise
                    error('MKMAP: Illegal value for title type (full/custom)')
            end
        end
        
        %% print image
        if printflag
            name=sprintf('%s_map%s_b%d_%s_%s_%s',CLUSTER,prjtag,boxx,printTypeTag,printtag,aexpn);
            %name=sprintf('%s/%s_map%s_b%d_%s.%s',printoutdir,CLUSTER,prjtag,boxx,printtag,'%s');
            
            printout_fig(gcf,name,'v')
        end
        
        
                
    end
end