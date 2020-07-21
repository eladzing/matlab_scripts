slice=mk_slice(s,ro,slind);
printag='s_zoom';

zmfac=0.05;

zmind=ceil((RVIR.*zmfac)./(boxx./0.7./256));

ll=zmind.*(1./0.7./256);
sid=[-ll ll];
sdind=129-zmind:128+zmind;


for projection = 1:3
  if plotproj(projection) 
      figure;
      
  if (clims(1)==clims(2))
     clims(1) = min(slice(:));
     clims(2) = max(slice(:));
  end
  
  %hold on
  imagesc(sid,sid,squeeze(log10(slice(sdind,sdind,projection))),clims);%axis equal tight;
%   resampled_v_x = imresize(squeeze(v_x(sdind,sdind,projection)), [diluted_len diluted_len], 'bilinear');
%   resampled_v_y = imresize(squeeze(v_y(sdind,sdind,projection)), [diluted_len diluted_len], 'bilinear');
%   h = streamslice(xxs,yys,squeeze(v_x(sdind,sdind,projection)),squeeze(v_y(sdind,sdind,projection)),1.5);
%   quiver(xxv,yyv,resampled_v_x,resampled_v_y, 0.9, 'w');
%   set(h,'color','k')
  zone_circles(zones);
  %hold off 
  
  switch projection
    case 1  % YZ
      xlabel('Z [Mpc/h]','Fontsize',14);
      ylabel('Y [Mpc/h]','Fontsize',14);
      prjtag='YZ'
    case 2  % ZX
      xlabel('X [Mpc/h]','Fontsize',14);
      ylabel('Z [Mpc/h]','Fontsize',14);
      prjtag='XZ'
    case 3  % XY
      xlabel('X [Mpc/h]','Fontsize',14);
      ylabel('Y [Mpc/h]','Fontsize',14);
      prjtag='XY'
  end  
  
  set(gca,'Ydir','normal');
  set(gca,'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1],'Fontsize',12)
  
  caxis(clims);
  colormap(avijet);
  bar=colorbar;
  set(get(bar,'Title'),'String','log(S/S_{vir})','Fontsize',14);


title(sprintf('%s Entropy, Thickness=%s Mpc/h',clustername,num2str(thick,3)),'Fontsize',14);


set(gcf,'Colormap',avijet);%'PaperOrientation','landscape')

if strcmp(pflag,'print')
   saveas(gcf,sprintf('%s/%s_map%s_b%d_%s.png',result_dir,clustername,prjtag,boxx,printag));
end

  end
end