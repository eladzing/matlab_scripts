function result = Vrcm_rvir(MPc)
% Utility function to load the temperature cube in spherical coordinates
% (See: cart2sphere)
%
% @param MPSec The required Mpc resolution (1, 2, 4, or 8)
% @returns  A 256^3 matrix of singles.
%

error(nargchk(1,1,nargin,'struct'));

global FILE_FORMAT_MAT
result = load(sprintf(FILE_FORMAT_MAT,'Vrcm_rvir',MPc),'res'); 
result = single(result.res);
