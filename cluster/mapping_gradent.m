figure;

slice=mk_slice(grads,ro,slind);
  
for projection = 1:3
  
  switch projection
    case 1  % YZ
      subplot(2,2,4);
    case 2  % ZX
      subplot(2,2,1);
    case 3  % XY
      subplot(2,2,3);
  end      
  
  if (clims(1)==clims(2))
     clims(1) = min(slice(:));
     clims(2) = max(slice(:));
  end
  
  hold on
  imagesc(side,side,squeeze(log10(slice(:,:,projection))),clims);%axis equal tight;
  resampled_v_x = imresize(squeeze(v_x(:,:,projection)), [diluted_len diluted_len], 'bilinear');
  resampled_v_y = imresize(squeeze(v_y(:,:,projection)), [diluted_len diluted_len], 'bilinear');
  h = streamslice(xxs,yys,squeeze(v_x(:,:,projection)),squeeze(v_y(:,:,projection)),1.5);
  quiver(xxv,yyv,resampled_v_x,resampled_v_y, 1.5, 'w');
  set(h,'color','k')
  zone_circles(zones);
  %hold off 
  
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

if strcmp(pflag,'print')
        saveas(gcf,sprintf('%s/%s_map3_b%d_%s.png',result_dir,clustername,boxx,printag));
end
