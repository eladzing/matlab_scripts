function mkmap_bare(cube,varargin)

% bare-bones map making

% arguments: 
% cube - data cube
% thick - half-thickness of slice in cells defualt=4
% velx,vely,velz - velocity field to plot - must already be density weighted
% dilute - vector field dilution ; defualt dilute=4 
% projection - which projection to plot  
% clims - limits of colorbar
% boxx - physical length of box
% weight - weighting cube for averaging
% circles - draw circles of given radii
% log - plot log10 of the data 

% defualt values 
weight=ones(size(cube));
thick=4;
dilute=4;
projection=3;
ncell=size(cube,1);  
boxx=1;
velflag=[0 0 0];
vx=0;vy=0;vz=0;
climflag=0;
logflag=0;

% read values
if nargin < 1
    error('ray_profile: not enough arguments')
elseif nargin > 23
    error('ray_profile: too many arguments')
end

for i=1:2:length(varargin)
    param=varargin{i};
    val=varargin{i+1};
    switch param
        case {'boxx','box'}
            boxx=val;
        case{'thick'}
            thick=val;
        case{'weightcube','wtcube','wtdatacube','wt','weight'}
            weight=val;
        case{'dilute'}
            dilute=val;
        case{'velx','vx'}
            vx=val;
            velflag(2)=1;
        case{'vely','vy'}
            vy=val;
            velflag(3)=1;
        case{'velz','vz'}
            vz=val;
            velflag(1)=1;
        case{'projection','proj'}
            proj=val;
            switch proj
                case {'yz','YZ','zy','ZY','x','X'}
                  projection=1;  
                case {'zx','ZX','xz','XZ','y','Y'}
                  projection=2;  
                case {'xy','XY','YX','yx','z','Z'}
                  projection=3;
                otherwise
                    error('mkmap_bare: illegal projection %s',proj);
            end
        case{'clims'}
            clims=val;
            climflag=1;
        case{'circles'}
            circles=val;
        case{'log','logdata'}
            logflag=val;
        otherwise
            error('mkmap_bare: illegal argument %s',param)
    end
end


slind=(floor(0.5*ncell)-thick+1):1:(ceil(0.5*ncell)+thick); 
zmind=0.5*ncell;

sidind=(0.5*ncell+1)-zmind:0.5*ncell+zmind;
actbox=zmind*(boxx/ncell);
side=[-actbox actbox];
zoomboxlen=length(sidind);

slice=mk_slice(cube,weight,slind);
clear cube weight `

if all(velflag)
    [v_x v_y]=mk_vfield(vx,vy,vz,ones(size(vx)),slind); 
    clear vx vy vz;
    diluted_len = length(1:dilute:zoomboxlen);
    diluted_jump = 2.*actbox/(diluted_len-1);
    notdiluted_jump = 2.*actbox/(zoomboxlen-1);
    [xxv yyv] = meshgrid(-actbox:diluted_jump:actbox, -actbox:diluted_jump:actbox);
    [xxs yys] = meshgrid(-actbox:notdiluted_jump:actbox,-actbox:notdiluted_jump:actbox);            
end

tickjump = actbox/4;
load('MyColormaps','avijet');

if ~climflag
    clims(1) = min(slice(:));
     clims(2) = max(slice(:));
end
hold on

%draw map
if logflag
    map=log10(squeeze(slice(:,:,projection)));
else
    map=squeeze(slice(:,:,projection));
end
imagesc(side,side,map,clims);%axis equal tight;
 
%draw velocity field
 if velflag 
      resampled_v_x = imresize(squeeze(v_x(:,:,projection)), [diluted_len diluted_len], 'bilinear');
      resampled_v_y = imresize(squeeze(v_y(:,:,projection)), [diluted_len diluted_len], 'bilinear');
      h = streamslice(xxs,yys,squeeze(v_x(sidind,sidind,projection)),squeeze(v_y(sidind,sidind,projection)),1.5);
      quiver(xxv,yyv,resampled_v_x,resampled_v_y, 1.3, 'w');
      set(h,'color','k')
  end
  
 % draw circles 
  if isempty(circles)>0
        draw_circles(circles);
  end
      
  hold off
  
  switch projection
    case 1  % YZ
      xlabelmine('$Z\, [\mathrm{Mpc/h}]$');
      ylabelmine('$Y\, [\mathrm{Mpc/h}]$');
    case 2  % ZX
      xlabelmine('$X\, [\mathrm{Mpc/h}]$');
      ylabelmine('$Z\, [\mathrm{Mpc/h}]$');
  
    case 3  % XY
      xlabelmine('$X\, [\mathrm{Mpc/h}]$');
      ylabelmine('$Y\, [\mathrm{Mpc/h}]$');
  
  end  
  
  set(gca,'Ydir','normal');
  set(gca,'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1],'XTick', -boxx/2:tickjump:boxx/2,'YTick', -boxx/2:tickjump:boxx/2, 'TickLength',[-0.015 -0.015],'Fontsize',12)
  set(gca,'XLim',[-boxx./2 boxx./2],'YLim',[-boxx./2 boxx./2])  
  
  caxis(clims);
  %colormap(avijet);
  %set(gcf,'Colormap',avijet);
end
  
  



