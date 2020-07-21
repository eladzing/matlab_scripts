function  plot_xrayProj_Chandra(proj) %,varargin )
%Plot_xrayProj_Chandra - - plot an x-ray projection as it would be observed
%  by Chandra

%units;

boxx=proj.box;
zemit=proj.zEmit;
%proj=yzA;   %xyA;   
cellsize=get_cellsize(boxx,'cm');

%inds=1:length(proj.ebins);

%%  get Chandra stuff
    fname='C:\\Users\\eladzing\\Documents\\cluster\\xray_code\\chandra.area';
    fid=fopen(fname);
    aa=fscanf(fid,'%g');
    area=(reshape(aa,[2 length(aa)/2]))';
    area2=interp1(area(:,1),area(:,2),proj.ebins);
    
    chanArea=zeros(size(proj.data));
    for i=1:size(chanArea,1)
        for j=1:size(chanArea,2)
            chanArea(i,j,:)=area2;
        end
    end

    
% numnerical factor: l^3/(4*pi*Dl^2)/(l/Da)^2 and conversion to arcsec^-2

sterad2arcsec=(180/pi*3600)^2; % convert steradian to arcsec^2
DaOverDl=(1+zemit).^-2;

fac=cellsize*DaOverDl^2/(4*pi*sterad2arcsec); 
%fac=cellsize*pi/(360*3600)^2*(1+zemit)^-4;
    
integrand=proj.data.*chanArea;
projI=trapz(proj.ebins,integrand,3).*fac;






figure
imagesc(boxx./2.*[-1 1],boxx./2.*[-1 1],log10(projI))
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
barTitle(hb,'$\mathrm{log(counts\,arcsec^{-2}\,s^{-1})}$')


end

