function mkmap(boxx,varargin)

% basic function for making a map of a cluster plus a velocity field.
% boxx -which boxx to map
% plotproj - which projection to plot, logical vector
% type - which value do we plot
%clims - limits of colorbar
%velflag - flag to create velocity field or not.
% vcm - center of mass velocity should be calculated in advance
% thick - the thickness of the slice in Mpc/h comoving
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
bartag='';

RVIR=get_rvir();
MVIR=get_mvir();
%VVIR=get_vvir();
TVIR=get_tvir();

% defualt values
%boxx=8;
plotproj=[1 1 1]; % plot all three projections
dilute=4;
marks=[];
vcm=[];
printoutdir=DEFAULT_PRINTOUT_DIR;
printtag='badtag';
printformat=[1 1];  % print eps and png
type='';
thick=-6   ;
normfactor=1;  % option to enter a normalization for the data
titletag='';
labeltype='full';
titletype='full';
rvin=logical([0 1 0]); % must be 3 component vector of 0 or 1 for plotting r500 rvir r200.
% default is to plot rvir alone.

%flags
cubeflag=false;
typeflag=false;
wtflag=false;
%labelflag=false;
titleflag=true;
printflag=false;
clim_flag=false; % no colorbar limits given, use min max
velflag=false;
%mvirboxflag=false;
logflag=false;
comoveflag=true;
vcubeflag(1:3)=false;
rogflag=false;
zoomFlag=false;
%% reading arguments

i=1;
while i<=length(varargin)
    switch varargin{i}
        case {'boxx','box'}   %which box to plot
            %           boxx=varargin{i};
        case{'data','datacube'} %plot a given cube;
            i=i+1;
            cubeflag=true;
            cube=varargin{i};
        case{'weightcube','wtcube','wtdatacube','wt'} % weight for cube
            i=i+1;
            wtflag=true;
            weight=varargin{i};
        case{'rocube','rog'} % rps cube essential for calculating contours
            i=i+1;
            rogflag=true;
            rog=varargin{i};
        case{'vxcube'}; % option to enter a pre-loaded velocity cube
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
        case{'datatype','type'} %prepare data cube of this type
            i=i+1;
            typeflag=true;
            type=varargin{i};
        case{'log','log10'}
            logflag=true;
        case{'nolog','linear'}
            logflag=false;
        case{'width','thick'} %thicknes of slice in Mpc/h comoving (like box)
            i=i+1;
            thick=varargin{i};
        case {'yz','YZ','zy','ZY'}
            plotproj=[1 0 0];
        case {'zx','ZX','xz','XZ'}
            plotproj=[0 1 0];
        case {'xy','XY','yx','YX'}
            plotproj=[0 0 1];
        case 'all'
            plotproj=[1 1 1];
        case {'proj','projection'}
            i=i+1;
            answer=varargin{i};
            if ischar(answer)
                switch answer
                    case {'yz','YZ','zy','ZY'}
                        plotproj=[1 0 0];
                    case {'zx','ZX','xz','XZ'}
                        plotproj=[0 1 0];
                    case {'xy','XY','yx','YX'}
                        plotproj=[0 0 1];
                    case 'all'
                        plotproj=[1 1 1];
                    otherwise
                        error('mkmap: unknown projection');
                end
            elseif isnumeric(answer)
                if length(answer)==1
                    switch answer
                        case 1
                            plotproj=[1 0 0];
                        case 2
                            plotproj=[0 1 0];
                        case 3
                            plotproj=[0 0 1];
                        otherwise
                            error('mkmap: unknown projection');
                    end
                elseif length(answer)==3
                    plotproj=answer;
                else
                    error('mkmap: unknown projection');
                end
            end
        case {'zoom','zoomBox'}
            i=i+1;
            zoomBox=varargin{i};
            zoomFlag=true;
        case {'clims','limits'}
            i=i+1;
            clim_flag=true;
            clims=varargin{i};
            if diff(clims)==0
                clim_flag=false;
            end
        case {'velflag','velocity','vfield'}
            velflag=true;
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
        case {'labels'}
            i=i+1;
            labeltype=varargin{i};
            %labelflag=(strcmp(l,'yes') || strcmp(answer,'on') || strcmp(answer));
        case {'title'}
            i=i+1;
            titletag=varargin{i};
            titleflag=~(strcmp(titletag,'no') || strcmp(titletag,'off') || strcmp(titletag,'none'));
        case {'titletype'}
            i=i+1;
            titletype=varargin{i};
        case {'bartag','bartitle'}
            i=i+1;
            bartag=varargin{i};
            %          case {'mvirbox'} % place a text of Mvir on the map
            %             mvirboxflag=true;
            %             bartag=varargin{i};
            %          case {'mvirbox'} % place a text of Mvir on the map
            %             mvirboxflag=true;
            %             strpos=varargin{i};  % position of tex            strpos=varargin{i};  % position of text in map coordinates
        case {'comoving','comove','comoving_coord'}
            comoveflag=true;
        case {'proper','proper_coord'}
            comoveflag=false;
        case {'print'}
            printflag=true;
        case {'noprint'}
            printflag=false;
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
            error('mkmap: illegal argument %s',varargin{i})
    end
    i=i+1;
end

if ~any(boxx==[1 2 4 8])
    error('mkmap: illegal box')
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



%bcson=sqrt(1.52e8.*T(boxx))./1e5; %speed of sound in km/sec

if(velflag || strcmpi(type,'flux'))
    if ~rogflag
        rog=RHOG(boxx);
        rogflag=true;
    end
    
    [hubX hubY hubZ] = hubble_flow(boxx,[0,0,0]);
    vx = Vx(boxx)+hubX-vcm(1);
    vy = Vy(boxx)+hubY-vcm(2);
    vz = Vz(boxx)+hubZ-vcm(3);
    vcubeflag(:)=true;
end




cm=[0,0,0];

if typeflag
    switch lower(type)
        case {'flux','massflux'}
            fluxnorm=(0.056.*(MVIR./1e13).^0.15*(1+zred)^2.25)./(4.*pi); % normalization for flux
            
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
            
            cube=rog.*vr.*rcube2./Mg200.*(1e5./3.0856e24.*3.15e7.*1e9)./fluxnorm ;
            weight=ones(size(cube));
            bartag='$\frac{\dot{M}_{gas}}{d\Omega}\,[\frac{\dot{M}_{gas}}{d\Omega}|_\mathrm{vir} ]$';
            %bartag='$\frac{\dot{M}/M_{gas}}{d\Omega} [\frac{\dot M}{M_{gas}}\left|_{\mathrm{virial}}]$';
            %bartag='$\dot{M}_{gas}$';
            
            printTypeTag='flux';
            
            clear vr rcube2 meshX meshY meshZ
            
        case {'density','gasdensity','rho','gas','rhogas'}
            if strcmpi(velflag,'vel')
                cube=rog;
            else
                cube=RHOG(boxx);
            end
            weight=ones(size(cube));
            logflag=true; %cube=log10(cube);
            bartag='$\log \rho_{gas}\,[\mathrm{M_\odot/Mpc^3}]$';
            printTypeTag='density';
        case {'dmdensity','rhodm','dm'}
            cube=RHODM(boxx);
            weight=ones(size(cube));
            logflag=true; %cube=log10(cube);
            bartag='$\log \rho_{DM}\,[\mathrm{M_\odot/Mpc^3}]$';
            printTypeTag='dmdensity';
        case {'totaldensity','rhotot','tot'}
            cube=RHOTOT(boxx);
            weight=ones(size(cube));
            logflag=true; %cube=log10(cube);
            bartag='$\log \rho_{tot}\,[\mathrm{M_\odot/Mpc^3}]$';
            printTypeTag='totdensity';
       
        case {'entropy','ent','k','s'}
            %SVIR=TVIR.*RVIR.^2./(3.*MVIR./(4.*pi)).^(2./3); %normalization for entropy
            cube=S(boxx).*f_ent; % in units of KeV cm^2 
            if strcmpi(velflag,'vel')
                weight=rog;
            else
                weight=RHOG(boxx);
            end
            logflag=true; %cube=log10(cube);
            bartag='$\log S\,[\mathrm{KeV\,cm^2}]$';
            printTypeTag='entropy';
        case {'temperature','temp','t'}
            cube=T(boxx);
            if strcmpi(velflag,'vel')
                weight=rog;
            else
                weight=RHOG(boxx);
            end
            logflag=true; %cube=log10(cube);
            bartag='$\log T\,[\mathrm{K}]$';
            printTypeTag='temperature';
        case {'pressure','press','p'}
            cube=T(boxx);
            if ~strcmpi(velflag,'vel')
                if ~rogflag
                    rog=RHOG(boxx);
                    rogflag=true;
                end
            end
            cube=cube.*rog./(TVIR*MVIR/(4*pi*RVIR^3/3));
            logflag=true; %cube=log10(cube);
            bartag='$\log P\,[\mathrm{P_{vir}}]$';
            printTypeTag='pressure';
        case {'metals','zia'}
            cube=ZIa(boxx);
            bartag='ZIa';
            printTypeTag='metalicityZIa';
        case {'iametals','zii'}
            cube=ZII(boxx);
            bartag='ZII';
            printTypeTag='metalicityZII';
        case {'metalicity','iimetals','ztot'}
            cube=ZIa(boxx)+ZII(boxx);
            bartag='Z';
            printTypeTag='metalicity';
        case {'potential','potent'}
            [pdm pgs , ~, ~]=poteng(boxx);
            cube=(pgs+pdm);
            %edot=ro.*pdot.*vol;
            clear pdm pgs ro
            bartag='$\phi\,[\mathrm{(km/sec)^2}]$';
            printTypeTag='potential_energy';
        case {'kinetic'}
            bartag='$E_{kin}\,[\mathrm{(km/sec)^2}]$';
            printTypeTag='kinetic';
        case {'poteng'}
            vol=(boxx./NCELL./hub).^3;
            
            [pdm pgs , ~, ro]=poteng(boxx);
            cube=(pgs+pdm).*ro.*vol;
            %edot=ro.*pdot.*vol;
            cube=-1*log10(abs(cube));
            clear pdm pgs pdot ro vol
            bartag='$E_{pot}\,[\mathrm{M_{\odot}(km/sec)^2}]$';
            printTypeTag='potential';
        case {'total_energy'}
        case {'cooling'}
            
        case {'velocity','vel'}
            if ~rogflag
                rog=RHOG(boxx);
                rogflag=true;
            end
            
            if any(~vcubeflag)
                [hubX hubY hubZ] = hubble_flow(boxx,[0,0,0]);
                vx = Vx(boxx)+hubX-vcm(1);
                vy = Vy(boxx)+hubY-vcm(2);
                vz = Vz(boxx)+hubZ-vcm(3);
                vcubeflag(:)=true;
            end
            
            cube=sqrt(vx.^2+vy.^2+vz.^2);
            weight=rog;
            
            %logflag=true; %cube=log10(cube);
            bartag='$v\,[\mathrm{km/sec}]$';
            printTypeTag='velocity';
        case {'vr','radialvel'}
            if ~rogflag
                rog=RHOG(boxx);
                rogflag=true;
            end
            
            if any(~vcubeflag)
                [hubX hubY hubZ] = hubble_flow(boxx,[0,0,0]);
                vx = Vx(boxx)+hubX-vcm(1);
                vy = Vy(boxx)+hubY-vcm(2);
                vz = Vz(boxx)+hubZ-vcm(3);
                vcubeflag(:)=true;
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
            %logflag=true; %cube=log10(cube);
            bartag='$v_r\,[\mathrm{km/sec}]$';
            printTypeTag='radvelocity';
        otherwise
            disp('unknown type');
            return
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
if(~comoveflag)
    side=side./(1+zred);
    actbox=actbox./(1+zred);
end
zoomBoxlen=length(sidind);

slice=mk_slice(cube,weight,slind);

if logflag
    slice=log10(slice);
end

% velocity field parameters (default filte=4)
if velflag
    [v_x v_y]=mk_vfield(vx,vy,vz,rog,slind);
    clear vx vy vz;
    diluted_len = length(1:dilute:zoomBoxlen);
    diluted_jump = 2.*actbox/(diluted_len-1);
    notdiluted_jump = 2.*actbox/(zoomBoxlen-1);
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
        if velflag
            resampled_v_x = imresize(squeeze(v_x(:,:,projection)), [diluted_len diluted_len], 'bilinear');
            resampled_v_y = imresize(squeeze(v_y(:,:,projection)), [diluted_len diluted_len], 'bilinear');
            h = streamslice(xxs,yys,squeeze(v_x(sidind,sidind,projection)),squeeze(v_y(sidind,sidind,projection)),1.5);
            quiver(xxv,yyv,resampled_v_x,resampled_v_y, 1.3, 'w');
            set(h,'color','k')
        end
        
        
        % make sure marked circles fit in box
        if ~isempty(marks)
            if(~comoveflag)
                markss=marks./(1+zred);
            else
                markss=marks;
            end
            marks_inbox=markss; %(marks<=(0.5.*boxx./hub));
            draw_circles(marks_inbox);
        end
        
        rv(1)=get_rvir(500);
        rv(2)=get_rvir();
        rv(3)=get_rvir(200);
        
        if(~comoveflag)
            rv=rv./(1+zred);
        end
        draw_circle_boxx(gcf,rv(rvin),'white');      % draw rvir
        
        
        hold off
        
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
        
        
        ticks=-boxx/(1+zred)/2:tickjump:boxx/(1+zred)/2;
        tickLabels=sprintf('%3.1f|',ticks);
        zl='|0.0|';
        k=strfind(tickLabels,zl);
        if k>0
            tickLabels=sprintf('%s%s%s',tickLabels(1:k-1),'|0|',tickLabels(k+length(zl):end));
        end
        
        set(gca,'Ydir','normal','Fontsize',12,...
            'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1],...
            'XTick',ticks ,'YTick',ticks,'XTickLabel',tickLabels,'YTickLabel',tickLabels,'TickLength',[-0.015 -0.015],...
            'XLim',side,'YLim',side)
        box on;
        
        if zoomFlag
            set(gca,'XLim',[-zoomBox zoomBox],'YLim',[-zoomBox zoomBox])
        end
        %else
        %    set(gca,'XLim',[-boxx./2 boxx./2],'YLim',[-boxx./2 boxx./2])
        %end
        
        
        %   if mvirboxflag
        %       ss=mk_mvir_string(MVIR);
        %       text(strpos(1),strpos(2),ss,'Interpreter','latex','fontsize',14);
        %   end
        
        
        caxis(clims);
        colormap(avijet);
        bar=colorbar;
        barTitle(bar,bartag);
        set(bar,'Fontsize',12)
        %set(get(bar,'Title'),'String',bartag,'Fontsize',12,'Interpreter','latex');
        set(gcf,'Colormap',avijet);
        %title(sprintf('%s %s, Thickness=%s Mpc/h',CLUSTER,type,num2str(thick,3)),'Fontsize',12,'Interpreter','latex');
        
        if titleflag
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
        
        %set(gcf,'Colormap',avijet);%'PaperOrientation','landscape')
        
        %printing:
        %Don't want more than 3 printing arugments altogether
        %numvarargs=length(varargin);
        %if numvarargs>3
        %    error('mkmap: too many printing arguments');
        %end
        
        %set defualt printing flag, printing tag, and output dir
        %defvals={'noprint', type ,'/home/eladzing/Ubuntu One/cluster/printout'};
        
        %assign the optional values
        %defvals(1:numvarargs)=varargin;
        
        % transfer to easy to use varaibles
        %[printflag printtag printoutdir]=defvals{:};
        
        
        if printflag
            name=sprintf('%s_map%s_b%d_%s_%s_%s',CLUSTER,prjtag,boxx,printTypeTag,printtag,aexpn);
            %name=sprintf('%s/%s_map%s_b%d_%s.%s',printoutdir,CLUSTER,prjtag,boxx,printtag,'%s');
            
            printout_fig(gcf,name,'dir',printoutdir,'v')
            
            
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
