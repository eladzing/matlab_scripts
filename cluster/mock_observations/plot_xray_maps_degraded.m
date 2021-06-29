%% plot sz maps and create text files. 

for k=[2 3 7 8 10 11 15 16] %1:length(xRayMap)
      %k=10; 
    new_env(xRayMap(k).cluster);
    
    global hub
    global NCELL
    global aexpn
    
    mapLen=8*NCELL;
    
    r500=get_rvir(500);
    r200=get_rvir(200);
    
    cl=get_cellsize(1,'mpc');
    
    cLen=ceil(min(2.1*r500/cl,mapLen/2));

    indx=mapLen/2-cLen+1:mapLen/2+cLen;
    xx=(indx-(mapLen/2+0.5)).*cl;

    projXY0=xRayMap(k).projXY(indx,indx);
    projYZ0=xRayMap(k).projYZ(indx,indx);
    projXZ0=xRayMap(k).projXZ(indx,indx);
    
   
    
    %% degrade 
    
    degScale = 500; 
    [yg,xg]=meshgrid(xx,xx);
    
    xxx=reshape(xg,[numel(xg) 1]);
    yyy=reshape(yg,[numel(yg) 1]);
   
    vvv=reshape(projXY0,[numel(projXY0) 1]);
    % pint=scatteredInterpolant(xxx,yyy,vvv);
    
    %[projXY, binsize, ~,~]= histogram2d(xxx,yyy,vvv,'len',degScale-1); 
    [projXY, binsize,~,~]= remap_fine2coarse(xxx,yyy,vvv,cl,'len',500,'wt',cl.^2);
    projXY=projXY(:,:,1)./prod(binsize);
    
    vvv=reshape(projYZ0,[numel(projYZ0) 1]);
    [projYZ, binsize,~,~]= remap_fine2coarse(xxx,yyy,vvv,cl,'len',500,'wt',cl.^2);
    projYZ=projYZ(:,:,1)./prod(binsize);
    
     
    vvv=reshape(projXZ0,[numel(projXZ0) 1]);
    [projXZ, binsize, xxlim,yylim]= remap_fine2coarse(xxx,yyy,vvv,cl,'len',500,'wt',cl.^2);
    projXZ=projXZ(:,:,1)./prod(binsize);
    
    
    xl=xxlim(1)+0.5.*binsize(1);
    %% plot images 
    global DEFAULT_PRINTOUT_DIR
    printoutdir=sprintf('%s/sz_data',DEFAULT_PRINTOUT_DIR);
    
    map = brewermap(256,'*RdBu');
    colormap(map);  
  
        
    
    figure
    imagesc(xxlim*hub,yylim*hub,log10(projXY));
    set(gca,'Ydir','Normal','Fontsize',14)
    bar=colorbar;set(bar,'Fontsize',12)
    barTitle(bar,'$\mathrm{log(counts\,arcsec^{-2}\,s^{-1})}$')
    caxis([-12 -2 ])
    draw_circle_boxx(gcf,r200,'white');      % draw rvir
    draw_circle_boxx(gcf,r500,'white');      % draw rvir
    draw_circle_boxx(gcf,2*r500,'white');      % draw rvir
    xlabelmine('$X\,[\mathrm{Mpc/h}]$');
    ylabelmine('$Y\,[\mathrm{Mpc/h}]$');
    titlemine(xRayMap(k).cluster);
    colormap(map);  
    axis square
    printout_fig(gcf,sprintf('%s_xrayMapXY_n%d_%s',xRayMap(k).cluster,degScale,aexpn),'subdir','sz_data','png')
    
    figure
    imagesc(xxlim*hub,yylim*hub,log10(projYZ));
    set(gca,'Ydir','Normal','Fontsize',14)
    bar=colorbar;set(bar,'Fontsize',12);barTitle(bar,'$\mathrm{log(counts\,arcsec^{-2}\,s^{-1})}$')
    draw_circle_boxx(gcf,r200,'white');      % draw rvir
    draw_circle_boxx(gcf,r500,'white');      % draw rvir
    draw_circle_boxx(gcf,2*r500,'white');      % draw rvir
    xlabelmine('$Z\,[\mathrm{Mpc/h}]$');
    ylabelmine('$Y\,[\mathrm{Mpc/h}]$');
    colormap(map);  titlemine(xRayMap(k).cluster); 
    caxis([-12 -2 ])
     axis square
    printout_fig(gcf,sprintf('%s_xrayMapZY_n%d_%s',xRayMap(k).cluster,degScale,aexpn),'subdir','sz_data','png')
    
    
    figure
    imagesc(xxlim*hub,yylim*hub,log10(projXZ));
    set(gca,'Ydir','Normal','Fontsize',14)
    bar=colorbar;set(bar,'Fontsize',12);barTitle(bar,'$\mathrm{log(counts\,arcsec^{-2}\,s^{-1})}$')
    draw_circle_boxx(gcf,r200,'white');      % draw rvir
    draw_circle_boxx(gcf,r500,'white');      % draw rvir
    draw_circle_boxx(gcf,2*r500,'white');      % draw rvir
    xlabelmine('$X\,[\mathrm{Mpc/h}]$');
    ylabelmine('$Z\,[\mathrm{Mpc/h}]$');
    colormap(map);   titlemine(xRayMap(k).cluster);
    caxis([-12 -2 ])
     axis square
    printout_fig(gcf,sprintf('%s_xrayMapXZ_n%d_%s',xRayMap(k).cluster,degScale,aexpn),'subdir','sz_data','png')
    
    %% write to file 
    
    
    head=sprintf('## bottom/left (Mpc): %s / %s , n_elements: %i x %i , cellSize: %s kpc, units: counts/sec/arcsec^2',...
        num2str(xl),num2str(xl),degScale,degScale,num2str(1000*binsize(1)));
    global HALO_PATH
    
    fnameXY=sprintf('%s/%s_xrayMap_xy_n%d_%s.dat',HALO_PATH,xRayMap(k).cluster,degScale,aexpn);
    fidXY=fopen(fnameXY,'w');
    fprintf(fidXY,'%s \n',head);
    fprintf(fidXY,'%12.5g \n',reshape(projXY',[numel(projXY) 1]));
    
    fnameYZ=sprintf('%s/%s_xrayMap_zy_n%d_%s.dat',HALO_PATH,xRayMap(k).cluster,degScale,aexpn);
    fidYZ=fopen(fnameYZ,'w');
    fprintf(fidYZ,'%s \n',head);
    fprintf(fidYZ,'%12.5g \n',reshape(projYZ',[numel(projYZ) 1]));
    
    fnameXZ=sprintf('%s/%s_xrayMap_xz_n%d_%s.dat',HALO_PATH,xRayMap(k).cluster,degScale,aexpn);
    fidXZ=fopen(fnameXZ,'w');
    fprintf(fidXZ,'%s \n',head);
    fprintf(fidXZ,'%12.5g \n',reshape(projXZ',[numel(projXZ) 1]));
    
    fclose('all');
    close all;
%     
%     
%     for i=1:length(indx)
%         for j=1:length(indx)
%             fprintf(fidXY,'%12.5g %12.5g %12.5g \n',xx(i),xx(j),projXY(i,j));
%             fprintf(fidYZ,'%12.5g %12.5g %12.5g \n',xx(i),xx(j),projYZ(i,j));
%             fprintf(fidXZ,'%12.5g %12.5g %12.5g \n',xx(i),xx(j),projXZ(i,j));
%             
%         end
%     end
  
 end