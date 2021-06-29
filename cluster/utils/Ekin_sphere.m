function [Ek Ekr] = Ekin_sphere(MPc)
% Utility function to load the kinetic energy cube cube in spherical coordinates
% (See: cart2sphere)
%
% @param MPSec The required Mpc resolution (1, 2, 4, or 8)
% @returns  A 256^3 matrix of singles.
%

error(nargchk(1,1,nargin,'struct'));

global FILE_FORMAT_SPHERE;
result = load(sprintf(FILE_FORMAT_SPHERE,'Ekin', MPc),'res'); 
result = single(result.res);
Ek=result(:,:,:,1); 
Ekr=result(:,:,:,2);
