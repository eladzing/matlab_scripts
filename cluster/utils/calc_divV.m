function div=calc_divV(Vxx,Vyy,Vzz,boxx)

% calculate the divergence of the cartesian velocity field
% assumes velocities are already 'fixed', i.e. with hubble flow
% and vcm corrected.

global NCELL
global hub

[meshX, meshY, meshZ] = meshgrid(1:NCELL, 1:NCELL, 1:NCELL);

%convert to center origin coordinates
        meshX = meshX - (NCELL+1)/2;% -cm(1);
        meshY = meshY - (NCELL+1)/2;% -cm(2);
        meshZ = meshZ - (NCELL+1)/2;% -cm(3);
        % Fix Units (to be in Mpc)
        meshX = meshX * ((boxx/hub)/NCELL);
        meshY = meshY * ((boxx/hub)/NCELL);
        meshZ = meshZ * ((boxx/hub)/NCELL);
    
        div=divergence(meshX,meshY,meshZ,Vxx,Vyy,Vzz);
        return