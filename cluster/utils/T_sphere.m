function result = T_sphere(MPc);
% Utility function to load the temperature cube in spherical coordinates
% (See: cart2sphere)
%
% @param MPSec The required Mpc resolution (1, 2, 4, or 8)
% @returns  A 256^3 matrix of singles.
%

error(nargchk(1,1,nargin,'struct'));

global FILE_FORMAT_SPHERE;
result = load(sprintf(FILE_FORMAT_SPHERE,'T', MPc),'res'); 
result = single(result.res);
