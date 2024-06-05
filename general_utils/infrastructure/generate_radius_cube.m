function rCube=generate_radius_cube(varargin)
%GENERATE_RADIUS_CUBE generate a grid with radial distances from the center
%of the cube stored at grid points 
%   Generate a 3d/2d grid containing the distance of each grid cell from
%   the center of cube. Input parmeter: 
%        nGrid - size of cube (can be an array of 1 2 or 3 values
%        leng - length of cube sides in physical units 
% 


cm=[0,0,0];
nGrid=256; 
leng=1;
cubeType='3d'; 
i=1;
while i<=length(varargin)
    switch lower(varargin{i})
        case {'2d','flat'}
            cubeType='2d';
            case {'3d'}
            cubeType='3d';
        case 'center'
            i=i+1;
            cm=varargin{i};
        case 'cube' % generate a cube with dimensions identical to a supplied cube 
            i=i+1;
            cube=varargin{i}; 
            nGrid=size(cube);
        case {'ng','ngrid','ncell','ncells'}
            i=i+1;
            nGrid=varargin{i};
        case {'length','size','side','sides'}
            i=i+1;
            leng=varargin{i};
        otherwise
            error('GENERATE_RADIUS_CUBE - Illegal argument: %s',varargin{i})
    end
    i=i+1;
end


%% asses the grid size
if all(size(nGrid)==[3 1])||all(size(nGrid)==[1 3])
    ng=nGrid;    % 3 different values given
elseif all(size(nGrid)==[2 1])||all(size(nGrid)==[1 2]) && strcmp(cubeType,'2d')
    ng=nGrid;     % 2 different values for a 2d grid 
elseif all(size(nGrid)==[1 1]) 
    % only one value for nGrid is given
    switch cubeType
        case '2d'
            ng=nGrid.*[1 1];
            case '3d'
            ng=nGrid.*[1 1 1];
    end
else
     error('GENERATE_RADIUS_CUBE - Illegal size argument: argument %s',nGrid);
end
            
%% assess the cube size in physical units.     
if all(size(leng)==[3 1])||all(size(leng)==[1 3])
    sz=leng;    % 3 different values given
elseif all(size(leng)==[2 1])||all(size(leng)==[1 2]) && strcmp(cubeType,'2d')
    sz=leng;     % 2 different values for a 2d grid 
elseif all(size(leng)==[1 1]) 
    % only one value for leng is given
    switch cubeType
        case '2d'
            sz=leng.*[1 1];
            case '3d'
            sz=leng.*[1 1 1];
    end
else
     error('GENERATE_RADIUS_CUBE - Illegal size argument: argument %s',leng);
end


[meshY, meshX, meshZ] = meshgrid(1:ng(1), 1:ng(2), 1:ng(3));

%find radial velocity component

%convert to center origin coordinates
meshX = meshX - (ng(1)+1)/2 -cm(1);
meshY = meshY - (ng(2)+1)/2 -cm(2);
meshZ = meshZ - (ng(3)+1)/2 -cm(3);
% Fix Units (to be in Mpc)
meshX = meshX .* sz(1)./ng(1);
meshY = meshY .* sz(2)./ng(2);
meshZ = meshZ .* sz(3)./ng(3);

rCube=sqrt(meshX.^2+meshY.^2+meshZ.^2) ; % r^2 cube in Mpc