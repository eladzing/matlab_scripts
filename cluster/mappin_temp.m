

global FILE_FORMAT;


cson=sqrt(1.52e8.*T(boxx))./1e5; %speed of sound in km/sec


ro=RHOG(boxx);
s=S(boxx);



if ~exist('cm')
    cm = [0,0,0];
end

vx = Vx(boxx);   vy = Vy(boxx);   vz = Vz(boxx);   

%%Vtot=sqrt(Vx.^2+Vy.^2+Vz.^2);
[meshY, meshX, meshZ] = meshgrid(1:size(vy,1), 1:size(vx,2), 1:size(vz,3));
%convert to center origin coordinates
meshX = meshX - (size(vx,1)+1)/2 -cm(1);
meshY = meshY - (size(vy,2)+1)/2 -cm(2);
meshZ = meshZ - (size(vz,3)+1)/2 -cm(3);
% Fix Units (to be in Mpc)
h=0.7;
meshX = meshX * ((boxx/h)/256);
meshY = meshY * ((boxx/h)/256);
meshZ = meshZ * ((boxx/h)/256);

%rcube2=meshX.^2+meshY.^2+meshZ.^2 ; % r^2 cube in Mpc
%vr=(vx.*meshX+vy.*meshY+vz.*meshZ)./sqrt(rcube2) ; %radial velocity 
%Mach=vr./cson; clear vr cson;


%% plot maps

dilute = 8;
blen=size(Mach,1);
thk=ceil(0.5.*thick./boxx.*blen);   %% thick is in Mpc
slind=(floor(0.5.*blen)-thk+1):1:(ceil(0.5.*blen)+thk);
side=[-boxx./2 boxx./2];
zone=[0.02 0.3 1.2];

[v_x v_y]=mk_vfield(vx,vy,vz,ro,slind) 
%clear vx vy vz
diluted_len = length(1:dilute:blen);
diluted_jump = boxx/(diluted_len-1);
notdiluted_jump = boxx/(blen-1);
[xxv yyv] = meshgrid(-boxx/2:diluted_jump:boxx/2, -boxx/2:diluted_jump:boxx/2);
[xxs yys] = meshgrid(-boxx/2:notdiluted_jump:boxx/2,-boxx/2:notdiluted_jump:boxx/2);
resampled_v_x = imresize(v_x, [diluted_len diluted_len], 'bilinear');
resampled_v_y = imresize(v_y, [diluted_len diluted_len], 'bilinear');
tickjump = boxx/8;
load_maxcolormap


%% plot entropy

figure;

slice=mk_slice(ss,ro,slind);
  
for projection = 1:3
  
  switch projection
    case 1  % YZ
      subplot(2,2,4);
    case 2  % ZX
      subplot(2,2,1);
    case 3  % XY
      subplot(2,2,3);
    end      
    
    imagesc(side,side,squeeze(log10(slice(:,:,projection))),clims);%axis equal tight;
    resampled_v_x = imresize(squeeze(v_x(:,:,projection)), [diluted_len diluted_len], 'bilinear');
    resampled_v_y = imresize(squeeze(v_y(:,:,projection)), [diluted_len diluted_len], 'bilinear');
    h = streamslice(xxs,yys,squeeze(v_x(:,:,projection)),squeeze(v_y(:,:,projection)),1.5);
    quiver(xxv,yyv,resampled_v_x,resampled_v_y, 1.5, 'w');
    set(h,'color','k')
    zone_circles(zone);
    
    switch projection
      case 1  % YZ
        xlabel('X [Mpc/h]','Fontsize',12);
      case 2  % ZX
        ylabel('Y [Mpc/h]','Fontsize',12);
      case 3  % XY
        xlabel('X [Mpc/h]','Fontsize',12);
        ylabel('Y [Mpc/h]','Fontsize',12);
      end      
      set(gca,'Ydir','normal');
    set(gca,'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1],'XTick', -boxx/2:tickjump:boxx/2,'YTick', -boxx/2:tickjump:boxx/2, 'TickLength',[-0.015 -0.015],'Fontsize',12)
  end
  
    
    
% XY
subplot(2,2,3);
imagesc(side,side,squeeze(log10(slice(:,:,3))),clims);%axis equal tight;
set(gca,'Ydir','normal');
xlabel('X [Mpc/h]','Fontsize',12);
ylabel('Y [Mpc/h]','Fontsize',12);
set(gca,'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1],'XTick', -boxx/2:tickjump:boxx/2,'YTick', -boxx/2:tickjump:boxx/2, 'TickLength',[-0.015 -0.015],'Fontsize',12)
zone_circles(zone);

resampled_v_x = imresize(squeeze(v_x(:,:,3)), [diluted_len diluted_len], 'bilinear');
resampled_v_y = imresize(squeeze(v_y(:,:,3)), [diluted_len diluted_len], 'bilinear');
h = streamslice(xxs,yys,squeeze(v_x(:,:,3)),squeeze(v_y(:,:,3)),1.5);
quiver(xxv,yyv,resampled_v_x,resampled_v_y, 1.5, 'w');
set(h,'color','k')






% ZX
subplot(2,2,1);
imagesc(side,side,squeeze(log10(slice(:,:,2))),clims);%axis equal tight;
set(gca,'Ydir','normal');
%%xlabel('X [Mpc/h]','Fontsize',12);
ylabel('Z [Mpc/h]','Fontsize',12);
set(gca,'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1],'XTick', -boxx/2:tickjump:boxx/2,'YTick', -boxx/2:tickjump:boxx/2, 'TickLength',[-0.015 -0.015],'Fontsize',12)
zone_circles(zone);

resampled_v_x = imresize(squeeze(v_x(:,:,2)), [diluted_len diluted_len], 'bilinear');
resampled_v_y = imresize(squeeze(v_y(:,:,2)), [diluted_len diluted_len], 'bilinear');
h = streamslice(xxs,yys,squeeze(v_x(:,:,2)),squeeze(v_y(:,:,2)),1.5);
quiver(xxv,yyv,resampled_v_x,resampled_v_y, 1.5, 'w');
set(h,'color','k')

% ZY
subplot(2,2,4);
imagesc(side,side,squeeze(log10(slice(:,:,1))),clims);%axis equal tight;
set(gca,'Ydir','normal');
xlabel('Z [Mpc/h]','Fontsize',12);
%ylabel('Z [Mpc/h]','Fontsize',12);grid;
set(gca,'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1],'XTick', -boxx/2:tickjump:boxx/2,'YTick', -boxx/2:tickjump:boxx/2, 'TickLength',[-0.015 -0.015],'Fontsize',12)
zone_circles(zone);

resampled_v_x = imresize(squeeze(v_x(:,:,1)), [diluted_len diluted_len], 'bilinear');
resampled_v_y = imresize(squeeze(v_y(:,:,1)), [diluted_len diluted_len], 'bilinear');
h = streamslice(xxs,yys,squeeze(v_x(:,:,1)),squeeze(v_y(:,:,1)),1.5);
quiver(xxv,yyv,resampled_v_x,resampled_v_y, 1.5, 'w');
set(h,'color','k')

% title box
subplot(2,2,2);
title(sprintf('%s Entropy Slicw Map Thickness=%s Mpc/h',clustername,num2str(thick,3)),'Fontsize',12);
	set(gca, 'XTick', [], 'YTick', []);
	%text(0.1,0.9, texttowrite, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    caxis(clims);
    colormap(maxcolormap);
    colorbar;


set(gcf,'Colormap',maxcolormap);%'PaperOrientation','landscape')



%%plot Mach
%figure; 
%    imagesc(side,side,(log10(abs(mean(Mach(:,:,slice),3))))');
%    %%load('MyColormaps','maxcolormap')
%    load_maxcolormap
%    set(gcf,'Colormap',maxcolormap)
%    bar=colorbar;
%    set(get(bar,'Title'),'String',sprintf('log(|v_r/c_s|)'),'Fontsize',12)
%    %caxis([-7 -2.5])

%    %set(gca,'Xdir','normal');
%    set(gca,'Ydir','normal');
    
%    xlabel('X [Mpc/h]','Fontsize',12);
%    ylabel('Y [Mpc/h]','Fontsize',12);grid;
%    title(sprintf( '%s Radial Velocity Mach Number (%d box)',clustername,num2str(boxx,2),'Fontsize',12));
%    set(gca,'fontsize',12)
      
    %if strcmp(pflag,'print')
%        saveas(gcf,sprintf('%s/%s_timebirds2_%s.png',result_dir,clustername,printag));
    %end

%end

