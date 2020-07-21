function [Ep Epdot] = Epot_sphere(MPc)
% Utility function to load the potential energy cube in spherical coordinates
% (See: cart2sphere)
%
% @param MPSec The required Mpc resolution (1, 2, 4, or 8)
% @returns  A 256^3 matrix of singles.
%

error(nargchk(1,1,nargin,'struct'));

global POT_FILE_FORMAT_SPHERE;
result = load(sprintf(POT_FILE_FORMAT_SPHERE,'Ep', MPc),'res'); 
Ep=result.res;
clear result
result = load(sprintf(POT_FILE_FORMAT_SPHERE,'Epdot', MPc),'res');
Epdot = result.res;

