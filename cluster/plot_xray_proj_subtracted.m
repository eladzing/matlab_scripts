global hub
global NCELL
units;
xl=[-0.6 -0.6+1.3];
yl=[-0.3 -0.3+1.3];

boxx=1;
inds=1:length(proj.ebins);
%proj=xyA; %yzA;   %xyA;
cellsize=boxx/hub/NCELL*Mpc;

projI=trapz(proj.ebins(inds),proj.data(:,:,inds),3).*cellsize;
%imaj2=10.^fittedmodel1.p2.*(rcube./get_rvir).^fittedmodel1.p1;
imaj2=(10.^fittedmodel.A).*(1+(rcube./get_rvir).^2.*fittedmodel.c.^2).^(0.5-3*fittedmodel.b);
figure

imagesc(boxx./2.*[-1 1],boxx./2.*[-1 1],log10(projI./imaj2))
%axis equal
set(gca,'Ydir','normal','Fontsize',14)
%xlim(xl);
%ylim(yl);

switch(lower(proj.projPlane))
    case{'xy','yx'}
        xlabelmine('$X\,[\mathrm{Mpc\,h^{-1}}]$')
        ylabelmine('$Y\,[\mathrm{Mpc\,h^{-1}}]$')
        
        
    case{'yz','zy'}
        xlabelmine('$Z\,[\mathrm{Mpc\,h^{-1}}]$')
        ylabelmine('$Y\,[\mathrm{Mpc\,h^{-1}}]$')
        
        
    case{'xz','zx'}
        xlabelmine('$X\,[\mathrm{Mpc\,h^{-1}}]$')
        ylabelmine('$Z\,[\mathrm{Mpc\,h^{-1}}]$')
end

caxis([-1 1])
cm=brewermap(256,'*RdBu');
colormap(cm)
hb=colorbar;
barTitle(hb,'log(Photon Flux/Profile)')