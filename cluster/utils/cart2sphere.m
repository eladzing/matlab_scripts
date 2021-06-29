function result = cart2sphere(cube,MAX_R,cm)
% Converts a cartesian cube to a spherical coordinates cube.
% Coordinates are (R,Phi,Theta) such that:
%    0.5   <= R     <= 128 (in jumps of 0.5)
%    -Pi/2 <= Phi   <= Pi/2 (in non-linear spacing to impose a constant ds
%                          per radius)
%    -Pi   <= Theta <= Pi
% See: sphere_grid, uni_sphere
%
% @param cube   The cartesian cube
% @param MAX_R  (optional) Specifies the maximum radius index (1-256)
%                           where 1 corresponds to 0.5 grid cell radius, and
%                               256 corresponds to 128 grid cell radius.
% @param cm     (optional) Center of coordinate system relative to the
%                          cube's center.
narginchk(1,3);

if ~exist('cm')
    cm = [0,0,0];
end

result = zeros(size(cube), 'single');

CUBE_SIZE = size(cube,1);
SPHERE_RES = CUBE_SIZE-1;
INTERP = 'linear';

[Sx, Sy, Sz] = uni_sphere(SPHERE_RES);
RR = [0.5:0.5:(CUBE_SIZE/2)];

if ~exist('MAX_R')
    MAX_R = length(RR);
end

for ridx = 1:MAX_R
    RSx = RR(ridx)*Sx+(128.5)+cm(1);
    RSy = RR(ridx)*Sy+(128.5)+cm(2);
    RSz = RR(ridx)*Sz+(128.5)+cm(3);
    
    result(ridx,:,:) = single(interp3(cube, RSx, RSy, RSz, INTERP));
end

end