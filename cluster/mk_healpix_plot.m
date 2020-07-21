%% construct healpix sampeling for cubes
global NCELL
global hub

if readFlag
    [meshX, meshY, meshZ] = meshgrid(1:NCELL, 1:NCELL, 1:NCELL);
    %convert to center origin coordinates
    meshX = meshX - (NCELL+1)/2;% -cm(1);
    meshY = meshY - (NCELL+1)/2;% -cm(2);
    meshZ = meshZ - (NCELL+1)/2;% -cm(3);
    
    cu(:,:,:,1)=T(1);
    cu(:,:,:,2)=T(2);
    cu(:,:,:,3)=T(4);
    cu(:,:,:,4)=T(8);
end
%% construct healpix on unit sphere
nHeal=64;
ss = HealpixGenerateSampling(nHeal, 'scoord');
rindx = HealpixGenerateSampling(nHeal, 'rindex');
%% create radial vector
dr=4/100;

r=dr:dr:4-dr;  %logspace(-1,log10(3.9),100); %  dr:dr:4-dr;
boxx=2.^(floor(log(2.*r)./log(2))+1);
boxx(boxx<1)=1;
ii=log(boxx)./log(2)+1;
arr=zeros(length(r),size(ss,1));

for i=1:length(r)
    C = SphToCartNew(ss,r(i)./hub);
    
    %% create X, Y ,Z positions
    
    
    % Fix Units (to be in Mpc)
    XX = meshX * ((boxx(i)/hub)/NCELL);
    YY = meshY * ((boxx(i)/hub)/NCELL);
    ZZ = meshZ * ((boxx(i)/hub)/NCELL);
    
    
    arr(i,:)=interp3(XX,YY,ZZ,squeeze(cu(:,:,:,ii(i))),C(:,1),C(:,2),C(:,3));
    
end


figure1=figure('Position',[1          49        1920         946]);
%        imagesc(1:lx,rr,log10(losMap1))
imagesc(1:size(arr,2),r./hub./get_rvir,log10(arr));
set(gca,'Ydir','Normal')
cmap=brewermap(256,'*Spectral');
colormap(cmap);
caxis([4 8])

bar=colorbar;
ylabelmine('$r/R_{\mathrm{vir}}$');
%barTitle(bar,'$\log(S)\,[\mathrm{K}]$')
set(gca,'Fontsize',14)
%xlabelmine('$\theta,\varphi$',16)
%     printout_fig(gcf,sprintf('cl%s_%s_losMap1',num2str(list(l)),aexpn))


ind2=switch_indexing(rindx);
figure1=figure('Position',[1          49        1920         946]);
%        imagesc(1:lx,rr,log10(losMap1))
imagesc(1:size(arr,2),r./hub./get_rvir,log10(arr(:,ind2)));
set(gca,'Ydir','Normal')
cmap=brewermap(256,'*Spectral');
colormap(cmap);
caxis([4 8])

bar=colorbar;
ylabelmine('$r/R_{\mathrm{vir}}$');
%barTitle(bar,'$\log(S)\,[\mathrm{K}]$')
set(gca,'Fontsize',14)
%xlabelmine('$\theta,\varphi$',16)
%     printout_fig(gcf,sprintf('cl%s_%s_losMap1',num2str(list(l)),aexpn))

