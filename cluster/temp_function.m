function slice=mk_slice(cube,weight,slind)
  
  cube=cube.*weight
  
  slice=zeros([size(cube,1),size(cube,2),3]);
  for projection = 1:3
    switch projection
      case 1
        slice(:,:,projection) = squeeze(sum(cube(slind,:,:),projection)./sum(weight(slind,:,:),projection));
      case 2
        slice(:,:,projection) = transpose(squeeze(sum(cube(:,slind,:),projection)./sum(weight(:,slind,:),projection)));
      case 3
        slice(:,:,projection) = transpose(squeeze(sum(cube(:,:,slind),projection)./sum(weight(:,:,slind),projection)));
    end
  end
  
end

function [v_x v_y]=mk_vfield(Vxx,Vyy,Vzz,weight,slind) 
  
  Vxx=Vxx.*weight;Vyy=Vyy.*weight;Vzz=Vzz.*weight;
  v_x=zeros([size(cube,1),size(cube,2),3])
  v_y=zeros([size(cube,1),size(cube,2),3])
  
  for projection = 1:3
        
    switch projection
      case 1
        v_x(:,:,projection) = squeeze(mean_nan(Vzz(slind,:,:),projection)./ mean_nan(weight(slind,:,:),projection));
        v_y(:,:,projection)  = squeeze(mean_nan(Vyy(slind,:,:),projection)./ mean_nan(weight(slind,:,:),projection));
            
      case 2
        v_y(:,:,projection)  = transpose(squeeze(mean_nan(Vzz(:,slind,:),projection)./ mean_nan(weight(:,slind,:),projection)));
        v_x(:,:,projection)  = transpose(squeeze(mean_nan(Vxx(:,slind,:),projection)./ mean_nan(weight(:,slind,:),projection)));
        
      case 3
        v_y(:,:,projection)  = transpose(squeeze(mean_nan(Vyy(:,:,slind),projection)./ mean_nan(weight(:,:,slind),projection)));
        v_x(:,:,projection)  = transpose(squeeze(mean_nan(Vxx(:,:,slind),projection)./ mean_nan(weight(:,:,slind),projection)));
        
      end
    end  
    
  end


if test==1 
%% plot gradient
figure;
printag='grad_s';

slice=mk_slice(ss,ones(size(s)),slind);
  
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
    case 1  % YZ;
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

% title box
subplot(2,2,2);
title(sprintf('%s Entropy Slicw Map Thickness=%s Mpc/h',clustername,num2str(thick,3)),'Fontsize',12);
set(gca, 'XTick', [], 'YTick', []);
%text(0.1,0.9, texttowrite, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
caxis(clims);
colormap(maxcolormap);
colorbar;


set(gcf,'Colormap',maxcolormap);%'PaperOrientation','landscape')

%if strcmp(pflag,'print')
%        saveas(gcf,sprintf('%s/%s_map3_b%d_%s.png',result_dir,clustername,boxx,printag));
    %end





%% plot Mach 
figure;
printag='mach';
slice=mk_slice(Mach,ones(size(Mach)),slind);
  
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

% title box
subplot(2,2,2);
title(sprintf('%s Entropy Slicw Map Thickness=%s Mpc/h',clustername,num2str(thick,3)),'Fontsize',12);
set(gca, 'XTick', [], 'YTick', []);
%text(0.1,0.9, texttowrite, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
caxis(clims);
colormap(maxcolormap);
colorbar;


set(gcf,'Colormap',maxcolormap);%'PaperOrientation','landscape')

%if strcmp(pflag,'print')
%        saveas(gcf,sprintf('%s/%s_map3_b%d_%s.png',result_dir,clustername,boxx,printag));
    %end

    
    
    %%plot flux
   figure;
printag='flux';
slice=mk_slice(flux,ones(size(flux)),slind);
  
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

% title box
subplot(2,2,2);
title(sprintf('%s Entropy Slicw Map Thickness=%s Mpc/h',clustername,num2str(thick,3)),'Fontsize',12);
set(gca, 'XTick', [], 'YTick', []);
%text(0.1,0.9, texttowrite, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
caxis(clims);
colormap(maxcolormap);
colorbar;


set(gcf,'Colormap',maxcolormap);%'PaperOrientation','landscape')

%if strcmp(pflag,'print')
%        saveas(gcf,sprintf('%s/%s_map3_b%d_%s.png',result_dir,clustername,boxx,printag));
    %end
end 




