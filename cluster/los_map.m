%% create radial los plot

list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
aa={'a1','a06'};
for l=15 ;%1:length(list)
    
    for n=1:2
        
        
        % read cube and convert to spherical coordinates
        boxx=8;
        new_env(list(l),aa{n});
        s=S(boxx);
        global NCELL
        global hub
        global aexpn
        ss=cart2sphere(s);
        
        % mask the directions
        sz=size(ss);
        mask=true([sz(2) sz(3)]);
        %rn=rand(size(mask));
        %mask(rn<0.5)=false;
        
        
        % extract radial los and restructure
        ly=sz(1);
        lx=sum(sum(mask));
        
        losMap1=zeros(ly,lx);
        losMap2=losMap1;
        for i=1:sz(1)
            rsh=reshape(ss(i,:,:),[1 lx]);
            losMap1(i,:)=rsh;
        end
        
        k=1;
        for i=1:sz(2)
            for j=1:sz(3)
                if ~mask(i,j)
                    continue
                end
                
                losMap2(:,k)=ss(:,i,j);
                k=k+1;
            end
        end
        
        rr=(1:256).*0.5.*boxx./hub./NCELL./get_rvir();
        
        
        figure1=figure('Position',[1          49        1920         946]);
        imagesc(1:lx,rr,log10(losMap1))
        set(gca,'Ydir','Normal')
        cmap=brewermap(256,'*Spectral');
        colormap(cmap);
        %caxis([4 8])
        bar=colorbar;
        ylabelmine('$r/R_{\mathrm{vir}}$');
        barTitle(bar,'$\log(T)\,[\mathrm{K}]$')
        set(gca,'Fontsize',14)
        xlabelmine('$\theta,\varphi$',16)
        printout_fig(gcf,sprintf('cl%s_%s_losMap1',num2str(list(l)),aexpn))
        
        figure2=figure('Position',[1          49        1920         946]);

        imagesc(1:lx,rr,log10(losMap2))
        set(gca,'Ydir','Normal')
        colormap(cmap);
        %caxis([4 8])
        bar=colorbar;
        set(bar,'Fontsize',16)
        ylabelmine('$r/R_{\mathrm{vir}}$',16);
        barTitle(bar,'$\log(T)\,[\mathrm{K}]$')
        set(gca,'Fontsize',16)
        xlabelmine('$\theta,\varphi$')
        printout_fig(gcf,sprintf('cl%s_%s_losMap2',num2str(list(l)),aexpn))
        
        %close all
        
    end
end