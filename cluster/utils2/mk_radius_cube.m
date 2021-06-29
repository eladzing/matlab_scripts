function rCube=mk_radius_cube(boxx,varargin)

global NCELL
global hub

hubl=hub;
cm=[0,0,0];
cubeFlag=false;
szFlag=false;
i=1;
sz=[1 1 1];
while i<=length(varargin)
    switch varargin{i}
        case 'cube'
            i=1+1;
            cube=varargin{i};
            cubeFlag=true;
        case 'center'
            i=i+1;
            cm=varargin{i};
        case {'ncell','size'}
            i=i+1;
            s0=varargin{i};
            if all(size(s0)==[3 1])||all(size(s0)==[1 3])
                sz=s0;
            elseif all(size(s0)==[1 1])
                sz=s0.*[1 1 1];
            else
                error('mk_radius_cube: Illegal size argument ')
            end
            
            szFlag=true;
        case {'hub','h'}
            i=i+1;
            hubl=varargin{i};
        otherwise
            error('mk_radius_cube: Illegal argument %s',varargin{i})
    end
    i=i+1;
end


if cubeFlag
    sz=[size(cube,1) size(cube,2) size(cube,3)];
elseif ~szFlag
    sz=NCELL.*[1 1 1];
end



[meshY, meshX, meshZ] = meshgrid(1:sz(1), 1:sz(2), 1:sz(3));

%find radial velocity component

%convert to center origin coordinates
meshX = meshX - (sz(1)+1)/2 -cm(1);
meshY = meshY - (sz(2)+1)/2 -cm(2);
meshZ = meshZ - (sz(3)+1)/2 -cm(3);
% Fix Units (to be in Mpc)
meshX = meshX * ((boxx/hubl)./sz(1));
meshY = meshY * ((boxx/hubl)./sz(2));
meshZ = meshZ * ((boxx/hubl)./sz(3));

rCube=sqrt(meshX.^2+meshY.^2+meshZ.^2) ; % r^2 cube in Mpc