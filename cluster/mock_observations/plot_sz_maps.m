%% plot sz maps and create text files. 

for k=10  %1:length(SZmap)
       
    new_env(SZmap(k).cluster);
    
    global hub
    global NCELL
    global aexpn
    
    mapLen=8*NCELL;
    
    r500=get_rvir(500);
    r200=get_rvir(200);
    
    cl=1/hub/NCELL;
    
    cLen=ceil(min(2.1*r500/cl,mapLen/2));

    indx=mapLen/2-cLen+1:mapLen/2+cLen;
    xx=(indx-(mapLen/2+0.5)).*cl;

    projXY=SZmap(k).projXY(indx,indx)./(cl.*1000).^2;
    projYZ=SZmap(k).projYZ(indx,indx)./(cl.*1000).^2;
    projXZ=SZmap(k).projXZ(indx,indx)./(cl.*1000).^2;
    
    
    %% plot images 
    global DEFAULT_PRINTOUT_DIR
    printoutdir=sprintf('%s/sz_data',DEFAULT_PRINTOUT_DIR);
    
    map = brewermap(256,'*RdBu');
    colormap(map);  
  
        
    
    figure
    imagesc(xx*hub,xx*hub,log10(projXY'));
    set(gca,'Ydir','Normal','Fontsize',14)
    bar=colorbar;set(bar,'Fontsize',12)
   barTitle(bar,'$\log(SZ)\,[\mathrm{kpc^{-2}}]$');
    caxis([-10 -5 ])
    draw_circle_boxx(gcf,r200,'white');      % draw rvir
    draw_circle_boxx(gcf,r500,'white');      % draw rvir
    draw_circle_boxx(gcf,2*r500,'white');      % draw rvir
    xlabelmine('$X\,[\mathrm{Mpc/h}]$');
    ylabelmine('$Y\,[\mathrm{Mpc/h}]$');
    titlemine(SZmap(k).cluster);
    colormap(map);  
    %printout_fig(gcf,sprintf('%s_szMapXY_%s',SZmap(k).cluster,aexpn),'subdir','sz_data')
    
    figure
    imagesc(xx*hub,xx*hub,log10(projYZ));
    set(gca,'Ydir','Normal','Fontsize',14)
    bar=colorbar;set(bar,'Fontsize',12)
    draw_circle_boxx(gcf,r200,'white');      % draw rvir
    draw_circle_boxx(gcf,r500,'white');      % draw rvir
    draw_circle_boxx(gcf,2*r500,'white');      % draw rvir
    xlabelmine('$Z\,[\mathrm{Mpc/h}]$');
    ylabelmine('$Y\,[\mathrm{Mpc/h}]$');
    colormap(map);  titlemine(SZmap(k).cluster); 
    barTitle(bar,'$\log(SZ)\,[\mathrm{kpc^{-2}}]$');
    caxis([-10 -5 ])
    %printout_fig(gcf,sprintf('%s_szMapZY_%s',SZmap(k).cluster,aexpn),'subdir','sz_data')
    
    
    figure
    imagesc(xx*hub,xx*hub,log10(projXZ'));
    set(gca,'Ydir','Normal','Fontsize',14)
     bar=colorbar;set(bar,'Fontsize',12)
    draw_circle_boxx(gcf,r200,'white');      % draw rvir
    draw_circle_boxx(gcf,r500,'white');      % draw rvir
    draw_circle_boxx(gcf,2*r500,'white');      % draw rvir
    xlabelmine('$X\,[\mathrm{Mpc/h}]$');
    ylabelmine('$Z\,[\mathrm{Mpc/h}]$');
    colormap(map);   titlemine(SZmap(k).cluster);
    barTitle(bar,'$\log(SZ)\,[\mathrm{kpc^{-2}}]$');
    caxis([-10 -5 ])
    printout_fig(gcf,sprintf('%s_szMapXZ_%s',SZmap(k).cluster,aexpn),'subdir','sz_data')
    
    %% write to file 
    
    tic 
    head=sprintf('## bottom/left (Mpc): %s / %s , n_elements: %i x %i , cellSize: %s kpc, units: kpc^-2',...
        num2str(xx(1)),num2str(xx(1)),length(indx),length(indx),num2str(1000*cl));
    global HALO_PATH
    
    fnameXY=sprintf('%s/%s_szMap_xy_%s.dat',HALO_PATH,SZmap(k).cluster,aexpn);
    fidXY=fopen(fnameXY,'w');
    fprintf(fidXY,'%s \n',head);
    fprintf(fidXY,'%12.5g \n',reshape(projXY,[numel(projXY) 1]));
    
    fnameYZ=sprintf('%s/%s_szMap_zy_%s.dat',HALO_PATH,SZmap(k).cluster,aexpn);
    fidYZ=fopen(fnameYZ,'w');
    fprintf(fidYZ,'%s \n',head);
    fprintf(fidYZ,'%12.5g \n',reshape(projYZ',[numel(projYZ) 1]));
    
    fnameXZ=sprintf('%s/%s_szMap_xz_%s.dat',HALO_PATH,SZmap(k).cluster,aexpn);
    fidXZ=fopen(fnameXZ,'w');
    fprintf(fidXZ,'%s \n',head);
    fprintf(fidXZ,'%12.5g \n',reshape(projXZ,[numel(projXZ) 1]));
    
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
  toc 
end