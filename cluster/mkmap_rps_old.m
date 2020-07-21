function mkmap_rps_old(boxx,clims,plotproj,thick,rps_type,rps_val,tag)  %velflag,dilute,marks,vcm,tag,varargin)

% basic function for making a map of a cluster plus a velocity field.  
% boxx -which boxx to map
% plotproj - which projection to plot, logical vector 
% rps_type - Which rps plot - by mass or for stripping radius? 
% clims - limits of colorbar
% rps_val - the value used to parameterize the calculation:
%            if rps_type is mvir then rps_val is mvir of the sat.
%            if rps_type is rstrip then rps_val is the .
%velflag - flag to create velocity field or not. 
% vcm - center of mass velocity should be calculated in advance 
% thick - the thickness of the slice in Mpc/h comoving
% marks - a vector of radii for drawing circles in units of Mpc
% dilute - vector field dilution a value <0 gives the defualt dilute=4

type='rps';
units;

switch rps_type
    case 'mvir'  
        %calculate the stripping radius for a given sat mv 
    [rv mv , ~, ~]=calculate_virials('mvir',rps_val);
    pval=(G/(8*pi)*(mv.*Ms)^2/(rv*pc*1e6)^4/Ms/km^2*Mpc^3);
    case 'rstrip'
    mm=1e12;
    [rv ~ , ~, ~]=calculate_virials('mvir',mm);
    alfa=rv./mm.^(1/3);
    pval=(G/(8*pi)/alfa^4/(rps_val^2))/Ms/km^2*M;
    otherwise
        error('ggg');
end

%bcson=sqrt(1.52e8.*T(boxx))./1e5; %speed of sound in km/sec


global CLUSTER;

%global HALO_PATH;
global NCELL;
global zred
global hub;
%load(sprintf('%s/virial%d_%s', HALO_PATH, 8,aexpn));

RVIR=get_rvir();
MVIR=get_mvir();
%VVIR=get_vvir();
TVIR=get_tvir();

[~ , ~, ~]=read_Mass_Profiles(RVIR);

% if (length(vcm)~=3)
%     vcm = VCM;
% end

weight=ones(NCELL,NCELL,NCELL);

%bcson=sqrt(1.52e8.*T(boxx))./1e5; %speed of sound in km/sec

% if(strcmpi(velflag,'vel')|| strcmpi(type,'flux'))
%     rog=RHOG(boxx);
%     [hubX hubY hubZ] = hubble_flow(boxx,[0,0,0]);
%     vx = Vx(boxx)+hubX-vcm(1);   
%     vy = Vy(boxx)+hubY-vcm(2);   
%     vz = Vz(boxx)+hubZ-vcm(3);   
%       resampled_v_x = imresize(squeeze(v_x(:,:,projection)), [diluted_len diluted_len], 'bilinear');
%       resampled_v_y = imresize(squeeze(v_y(:,:,projection)), [diluted_len diluted_len], 'bilinear');
%       h = streamslice(xxs,yys,squeeze(v_x(sidind,sidind,projection)),squeeze(v_y(sidind,sidind,projection)),1.5);
% end

cm=[0,0,0];


switch lower(type)
    
    case {'rampressure','rps','ram'} 
        rog=RHOG(boxx); 
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
        fac=G*Ms/(pc*1e6)/km^2;
        cube=rog.*mtot./rcube.*fac;
        %cube=1./sqrt((rog.*mtot./rcube.*fac)./pval);
        clear rcube mtot
        
        weight=ones(size(cube));
        cube=log10(cube);
        bartag='$\log \rho_{tot}v^2\,[\mathrm{km^2/sec^2/Mpc^3}]$';        
    otherwise
        disp('unknown type');
        return
end
        
%prepare map parameters 

%define slice index 
thk=ceil(0.5.*thick./boxx.*NCELL);   %% thick is in comoving Mpc
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
 cside=side(1):diff(side)/(NCELL-1):side(2);
 zoomboxlen=length(sidind);

slice=mk_slice(cube,weight,slind);
  
% % velocity field parameters 
% if dilute<=0 
%     dilute = 4;
% end
% if strcmpi(velflag,'vel')
%     [v_x v_y]=mk_vfield(vx,vy,vz,rog,slind); 
%     clear vx vy vz;
%     diluted_len = length(1:dilute:zoomboxlen);
%     diluted_jump = 2.*actbox/(diluted_len-1);
%     notdiluted_jump = 2.*actbox/(zoomboxlen-1);
%     [xxv yyv] = meshgrid(-actbox:diluted_jump:actbox, -actbox:diluted_jump:actbox);
%     [xxs yys] = meshgrid(-actbox:notdiluted_jump:actbox,-actbox:notdiluted_jump:actbox);            
% end
% diluted_len], 'bilinear');
%       h = streamslice(xxs,yys,squeeze(v_x(sidind,sidind,projection)),squeeze(v_y(sidind,sidind,projection)),1.5);
%       quiver(xxv,yyv,resampled_v_x,resampled_v_y, 1.3, 'w');
%       set(h,'color','k')
%   end
   
    
%   % make sure marked circles fit in box
%   if ~isempty(marks)
%         marks_inbox=marks; %(marks<=(0.5.*boxx./hub));

tickjump = actbox/4;
load('MyColormaps','avijet');



% draw map
for projection = 1:3
  if plotproj(projection)  
      figure;
  
  
      
  if (clims(1)==clims(2))
     clims(1) = min(slice(:));
     clims(2) = max(slice(:));
  end
  
  hold on
  islice=squeeze(slice(:,:,projection));
  imagesc(side,side,islice,clims);%axis equal tight;
  switch rps_type
      case 'mvir'
        rstr=1./sqrt((10.^islice)./pval);
        conts=[0.01 0.05 0.1 0.2 0.5 0.75 1.0];
      case 'rstrip'
        rstr=log10(((10.^islice)./pval).^(3/2));  
        conts=[9 10 11 12 13];
      otherwise
          error('kkkk');
  end
        C=contour(cside,cside,rstr,conts);
        clabel(C);
%   if (strcmpi(velflag,'vel')) 
%       resampled_v_x = imresize(squeeze(v_x(:,:,projection)), [diluted_len diluted_len], 'bilinear');
%       resampled_v_y = imresize(squeeze(v_y(:,:,projection)), [diluted_len diluted_len], 'bilinear');
%       h = streamslice(xxs,yys,squeeze(v_x(sidind,sidind,projection)),squeeze(v_y(sidind,sidind,projection)),1.5);
%       quiver(xxv,yyv,resampled_v_x,resampled_v_y, 1.3, 'w');
%       set(h,'color','k')
%   end
   
    
%   % make sure marked circles fit in box
%   if ~isempty(marks)conts=[0.01 0.05 0.1 0.2 0.5 0.75 1.0];
%         marks_inbox=marks; %(marks<=(0.5.*boxx./hub));
%         draw_circles(marks_inbox);
%   end
  draw_circle_boxx(gcf,get_rvir(),'white');      % draw rvir 

    
  hold off 
  
  switch projection
    case 1  % YZ
      xlabelmine('$Z\, [\mathrm{Mpc/h}]$');
      ylabelmine('$Y\, [\mathrm{Mpc/h}]$');
      prjtag='YZ';
    case 2  % ZX
      xlabelmine('$X\, [\mathrm{Mpc/h}]$');
      ylabelmine('$Z\, [\mathrm{Mpc/h}]$');
      prjtag='XZ';
    case 3  % XY
      xlabelmine('$X\, [\mathrm{Mpc/h}]$');
      ylabelmine('$Y\, [\mathrm{Mpc/h}]$');
      prjtag='XY';
   end  
  
    
  set(gca,'Ydir','normal');
  set(gca,'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1],'XTick', -boxx/2:tickjump:boxx/2,'YTick', -boxx/2:tickjump:boxx/2, 'TickLength',[-0.015 -0.015],'Fontsize',12)
  set(gca,'XLim',[-boxx./2 boxx./2],'YLim',[-boxx./2 boxx./2])  
  
  caxis(clims);
  colormap(avijet);
  bar=colorbar;
  set(get(bar,'Title'),'String',bartag,'Fontsize',12,'Interpreter','latex');
  set(gcf,'Colormap',avijet);
  %title(sprintf('%s %s, Thickness=%s Mpc/h',CLUSTER,type,num2str(thick,3)),'Fontsize',12,'Interpreter','latex');
  titlemine(sprintf('%s %s, %s, $z=%3.2g$',CLUSTER,type,tag,zred));

%set(gcf,'Colormap',avijet);%'PaperOrientation','landscape')

% %printing:
% %Don't want more than 3 printing arugments altogether
% numvarargs=length(varargin);
% if numvarargs>3
%     error('mkmap: too many printing arguments');
% end
% 
% %set defualt printing flag, printing tag, and output dir
% defvals={'noprint', type ,'/home/eladzing/Ubuntu One/cluster/printout'};
% 
% %assign the optional values
% defvals(1:numvarargs)=varargin;
% 
% % transfer to easy to use varaibles
% [pflag ptag printoutdir]=defvals{:};
% 
% 
% if strcmp(pflag,'print')
%     name=sprintf('%s/%s_map%s_b%d_%s.%s',printoutdir,CLUSTER,prjtag,boxx,ptag,'%s');
%     exportfig(gcf,sprintf(name,'png'),'format','png');
%     exportfig(gcf,sprintf(name,'eps'));
%     %saveas(gcf,sprintf('%s/%s_map%s_b%d_%s.png',printoutdir,CLUSTER,prjtag,boxx,ptag));    
% end

  end
end
end