function rcube=mk_rcube(boxx,cube,cm) %% returns cube of radial positions
global hub;
global aexp

if ~exist('cm')
    cm = [0,0,0];
end

if max(size(cube))==1;
    lx=cube;
    ly=cube;
    lz=cube;
else
    lx=size(cube,1);
    ly=size(cube,2);
    lz=size(cube,3);
end

[meshY, meshX, meshZ] = meshgrid(1:lx, 1:ly, 1:lz);

%convert to center origin coordinates
meshX = meshX - (lx+1)/2 -cm(1);
meshY = meshY - (ly+1)/2 -cm(2);
meshZ = meshZ - (lz+1)/2 -cm(3);
% Fix Units (to be in Mpc)
h=hub;
meshX = meshX * ((aexp*boxx/h)/256);
meshY = meshY * ((aexp*boxx/h)/256);
meshZ = meshZ * ((aexp*boxx/h)/256);

rcube=sqrt(meshX.^2+meshY.^2+meshZ.^2); % r cube in Mpc