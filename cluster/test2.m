cm=[0 0 0];
boxx=8;
hub=0.7;
NCELL=256;

[meshY, meshX, meshZ] = meshgrid(1:256, 1:256, 1:256);
        %convert to center origin coordinates
        meshX = meshX - (256+1)/2 -cm(1);
        meshY = meshY - (256+1)/2 -cm(2);
        meshZ = meshZ - (256+1)/2 -cm(3);
        % Fix Units (to be in Mpc)
        meshX = meshX * ((boxx/hub)/NCELL);
        meshY = meshY * ((boxx/hub)/NCELL);
        meshZ = meshZ * ((boxx/hub)/NCELL);

        rcube=sqrt(meshX.^2+meshY.^2+meshZ.^2); % r cube in Mpc
      