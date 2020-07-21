%% plot vr vs r for all the gas in the cluster

list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

epoch='a06';
len=512;

bigbird=zeros(len,len,length(list));


for j=1:length(list)
    new_env(list(j),epoch)
    
    global hub
    global NCELL
    global aexpn
    global CLUSTER
    
    mask1=true(NCELL,NCELL,NCELL);
    mask2=mask1;
    ll=floor(NCELL/4)+1:floor(3*NCELL/4);
    mask2(ll,ll,ll)=false;
    
    cm=[0 0 0];
    for i=1:4
        boxx=2^(i-1);
        rhog=RHOG(boxx);
        [vx, vy, vz] = get_velocities(boxx);
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
        
        rcube=sqrt(meshX.^2+meshY.^2+meshZ.^2) ; % r^2 cube in Mpc
        vr=(vx.*meshX+vy.*meshY+vz.*meshZ)./rcube ; %radial velocity
        vol=(boxx/hub/NCELL)^3;
        switch boxx
            case 1
                r1=rcube(mask1);
                ro1=rhog(mask1).*vol;
                v1=vr(mask1);
            case 2
                r2=rcube(mask2);
                ro2=rhog(mask2).*vol;
                v2=vr(mask2);
            case 4
                r4=rcube(mask2);
                ro4=rhog(mask2).*vol;
                v4=vr(mask2);
            case 8
                r8=rcube(mask2);
                ro8=rhog(mask2).*vol;
                v8=vr(mask2);
            otherwise
                error('Whaaaaat?')
        end
    end
    
    vr=cat(1,v1,v2,v4,v8)./get_vvir;
    ro=cat(1,ro1,ro2,ro4,ro8);
    rr=cat(1,r1,r2,r4,r8)./get_rvir;
    
    % [bird, ~, xxlim,yylim]= histogram2d(rr,vr,ro,'len',512);
    %
    % figure
    % imagesc(xxlim,yylim,log10(squeeze(bird(:,:,1))./sum(ro)));
    % set(gca,'Ydir','Normal')
    % cmap=flipud(linear_bgyw_15_100_c68_n256);% viridis(256);      %  brewermap(256,'Reds');
    % colormap(cmap);
    % hb=colorbar;
    %
    % caxis([-6 -3.5])
    %
    % barTitle(hb,'$\log(M/M_{\mathrm{gas}})$');
    % set(hb,'Fontsize',14)
    %
    % set(gca,'Fontsize',14)
    % hold on
    % plot(xxlim,get_vvir.*[1 1],'--k','linewidth',1.5,'DisplayName','$V_{\mathrm{vir}}$');
    % plot(xxlim,-get_vvir.*[1 1],'--k','linewidth',1.5,'DisplayName','$V_{\mathrm{vir}}$');
    % xlabelmine('$r/R_{\mathrm{vir}}$');
    % ylabelmine('$v_r\,[\mathrm{km\,sec^{-1}}]$');
    
    
    %% create bird
    [bird, ~, xxlim,yylim]= histogram2d(log10(rr),vr,ro,'len',len-1,'xxlim',log10([0.05 5.5]),'yylim',[-1.75 1.75]);
    
    bigbird(:,:,j)=squeeze(bird(:,:,1))./sum(ro);
    
    %% plot
    hf=figure;
    imagesc(xxlim,yylim,log10(squeeze(bird(:,:,1))./sum(ro)));
    set(gca,'Ydir','Normal')
    cmap=flipud(linear_bgyw_15_100_c68_n256);% viridis(256);      %  brewermap(256,'Reds');
    colormap(cmap);
    hb=colorbar;
    grid
    caxis([-6 -3.5])
    
    barTitle(hb,'$\log(M/M_{\mathrm{gas}})$');
    set(hb,'Fontsize',14)
    
    set(gca,'Fontsize',14)
    %hold on
    %plot(xxlim,get_vvir.*[1 1],'--k','linewidth',1.5,'DisplayName','$V_{\mathrm{vir}}$');
    %plot(xxlim,-get_vvir.*[1 1],'--k','linewidth',1.5,'DisplayName','$V_{\mathrm{vir}}$');
    xlabelmine('$\log(r/R_{\mathrm{vir}})$');
    ylabelmine('$v_r/V_{\mathrm{vir}}$');
    titlemine(CLUSTER)
    
    fname=sprintf('%s_r_vr_bird_%s.pdf',CLUSTER,aexpn);
    printout_fig(gcf,fname,'nofig')
end