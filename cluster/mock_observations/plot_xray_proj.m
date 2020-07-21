global hub
global NCELL
units;
%xl=[-0.6 -0.6+1.3];
%yl=[-0.3 -0.3+1.3];

boxx=proj.box;

%proj=yzA;   %xyA;

cellsize=boxx/hub/NCELL*Mpc;
inds=1:length(proj.ebins);
projI=trapz(proj.ebins(inds),proj.data(:,:,inds),3).*cellsize;
figure
imagesc(boxx./2.*[-1 1],boxx./2.*[-1 1],log(projI))
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


cm=brewermap(256,'*Greys');
colormap(cm)
hb=colorbar;
barTitle(hb,'$\mathrm{log(counts\,cm^{-2}\,s^{-1})}$')