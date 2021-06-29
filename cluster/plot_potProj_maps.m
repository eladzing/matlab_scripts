%% plot sz maps and create text files. 

for k=1:length(potMap)
       
    new_env(potMap(k).cluster);
    
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

    projXY0=potMap(k).projXY(indx,indx);%./(cl.*1000);
    projYZ0=potMap(k).projYZ(indx,indx);%./(cl.*1000);
    projXZ0=potMap(k).projXZ(indx,indx);%./(cl.*1000);
    
     %% degrade 
    
    degScale = 500; 
    [yg,xg]=meshgrid(xx,xx);
    
    xxx=reshape(xg,[numel(xg) 1]);
    yyy=reshape(yg,[numel(yg) 1]);
   
    vvv=reshape(projXY0,[numel(projXY0) 1]);
    [projXY, binsize, xxlim,yylim]= histogram2d(xxx,yyy,vvv,'len',degScale-1); 
    projXY=projXY(:,:,1)'./prod(binsize.*1000);
    
    vvv=reshape(projYZ0,[numel(projYZ0) 1]);
    [projXZ, binsize, xxlim,yylim]= histogram2d(xxx,yyy,vvv,'len',degScale-1); 
    projXZ=projXZ(:,:,1)'./prod(binsize.*1000);
    
    vvv=reshape(projXZ0,[numel(projXZ0) 1]);
    [projYZ, binsize, xxlim,yylim]= histogram2d(xxx,yyy,vvv,'len',degScale-1); 
    projYZ=projYZ(:,:,1)'./prod(binsize.*1000);
    
    %% plot images 
    global DEFAULT_PRINTOUT_DIR
    printoutdir=sprintf('%s/sz_data',DEFAULT_PRINTOUT_DIR);
    
    map = brewermap(256,'*RdBu');
    colormap(map);  
  
        
    bartag='$\log(|\Phi|)\,[\mathrm{((km/sec)^2\,kpc^{-2}}]$';
    figure
    imagesc(xx*hub,xx*hub,log10(abs(projXY0')));
    set(gca,'Ydir','Normal','Fontsize',14)
    bar=colorbar;set(bar,'Fontsize',12)
    barTitle(bar,bartag);
    caxis([-10 -5 ])
    draw_circle_boxx(gcf,r200,'white');      % draw rvir
    draw_circle_boxx(gcf,r500,'white');      % draw rvir
    draw_circle_boxx(gcf,2*r500,'white');      % draw rvir
    xlabelmine('$X\,[\mathrm{Mpc/h}]$');
    ylabelmine('$Y\,[\mathrm{Mpc/h}]$');
    titlemine(potMap(k).cluster);
    colormap(map);  
    printout_fig(gcf,sprintf('%s_potMapXY_%s',potMap(k).cluster,aexpn),'subdir','sz_data')
    
    figure
    imagesc(xx*hub,xx*hub,log10(abs(projYZ)));
    set(gca,'Ydir','Normal','Fontsize',14)
    bar=colorbar;set(bar,'Fontsize',12)
    draw_circle_boxx(gcf,r200,'white');      % draw rvir
    draw_circle_boxx(gcf,r500,'white');      % draw rvir
    draw_circle_boxx(gcf,2*r500,'white');      % draw rvir
    xlabelmine('$Z\,[\mathrm{Mpc/h}]$');
    ylabelmine('$Y\,[\mathrm{Mpc/h}]$');
    colormap(map);  titlemine(potMap(k).cluster); 
    barTitle(bar,bartag);
    caxis([-10 -5 ])
    printout_fig(gcf,sprintf('%s_potMapZY_%s',potMap(k).cluster,aexpn),'subdir','sz_data')
    
    
    figure
    imagesc(xx*hub,xx*hub,log10(abs(projXZ')));
    set(gca,'Ydir','Normal','Fontsize',14)
     bar=colorbar;set(bar,'Fontsize',12)
    draw_circle_boxx(gcf,r200,'white');      % draw rvir
    draw_circle_boxx(gcf,r500,'white');      % draw rvir
    draw_circle_boxx(gcf,2*r500,'white');      % draw rvir
    xlabelmine('$X\,[\mathrm{Mpc/h}]$');
    ylabelmine('$Z\,[\mathrm{Mpc/h}]$');
    colormap(map);   titlemine(potMap(k).cluster);
    barTitle(bar,bartag);
    caxis([-10 -5 ])
    printout_fig(gcf,sprintf('%s_potMapXZ_%s',potMap(k).cluster,aexpn),'subdir','sz_data')
    
    %% write to file 
    
    tic 
    head=sprintf('## bottom/left (Mpc): %s / %s , n_elements: %i x %i , cellSize: %s kpc, units: (km/sec)^2 kpc^-2',...
        num2str(xx(1)),num2str(xx(1)),length(indx),length(indx),num2str(1000*cl));
    global HALO_PATH
    
    fnameXY=sprintf('%s/%s_potMap_xy_%s.dat',HALO_PATH,potMap(k).cluster,aexpn);
    fidXY=fopen(fnameXY,'w');
    fprintf(fidXY,'%s \n',head);
    fprintf(fidXY,'%12.5g \n',reshape(projXY,[numel(projXY) 1]));
    
    fnameYZ=sprintf('%s/%s_potMap_zy_%s.dat',HALO_PATH,potMap(k).cluster,aexpn);
    fidYZ=fopen(fnameYZ,'w');
    fprintf(fidYZ,'%s \n',head);
    fprintf(fidYZ,'%12.5g \n',reshape(projYZ',[numel(projYZ) 1]));
    
    fnameXZ=sprintf('%s/%s_potMap_xz_%s.dat',HALO_PATH,potMap(k).cluster,aexpn);
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