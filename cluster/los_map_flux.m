%% create radial los plot

list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
aa={'a1','a06'};
for l=15 ;%1:length(list)
    
    for n=2 %1;%:2
        
        
        % read cube and convert to spherical coordinates
        %boxx=8;
        new_env(list(l),aa{n});
        global NCELL
        global hub
        global aexpn

        s1=new_flux(1);
        s2=new_flux(2);
        s4=new_flux(4);
        s8=new_flux(8);

        ss=cat(1,s1(10:end-1,:,:),s2(129:end-1,:,:),s4(129:end-1,:,:),s8(129:end-1,:,:));
        r1=(10:255).*1;
        r2=(129:255).*2;
        r4=(129:255).*4;
        r8=(129:255).*8;
        rr=cat(2,r1,r2,r4,r8)./hub./NCELL./get_rvir();
        %%ss=cart2sphere(s);
        
        
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
        
        %rr=(1:256).*boxx./hub./NCELL./get_rvir();
        cmap=brewermap(256,'*Spectral');
        
        figure1=figure('Position',[1          49        1920         946]);
        imagesc(1:lx,rr,(losMap1))
        set(gca,'Ydir','Normal')
        colormap(cmap);
        caxis([-0.1 0.1])
        bar=colorbar;
        ylabelmine('$r/R_{\mathrm{vir}}$');
        barTitle(bar,'$\log(T)\,[\mathrm{K}]$')
        set(gca,'Fontsize',14)
        xlabelmine('$\theta,\varphi$',16)
        %printout_fig(gcf,sprintf('cl%s_%s_losMap1',num2str(list(l)),aexpn))
        
        figure2=figure('Position',[1          49        1920         946]);

        imagesc(1:lx,rr,(losMap2))
        set(gca,'Ydir','Normal')
        colormap(cmap);
        caxis([-0.1 0.1])
        bar=colorbar;
        set(bar,'Fontsize',16)
        ylabelmine('$r/R_{\mathrm{vir}}$',16);
        barTitle(bar,'$\log(T)\,[\mathrm{K}]$')
        set(gca,'Fontsize',16)
        xlabelmine('$\theta,\varphi$')
        %printout_fig(gcf,sprintf('cl%s_%s_losMap2',num2str(list(l)),aexpn))
        
        %close all
        
    end
end