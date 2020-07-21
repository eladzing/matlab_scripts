function [Potgs Potdm Potdot] = Pot_sphere(MPc)
% Utility function to load the potential energy cube in spherical coordinates
% (See: cart2sphere)
%
% @param MPSec The required Mpc resolution (1, 2, 4, or 8)
% @returns  A 256^3 matrix of singles.
%

error(nargchk(1,1,nargin,'struct'));

global POT_FILE_FORMAT_SPHERE;
result = load(sprintf(POT_FILE_FORMAT_SPHERE,'Potgs', MPc),'res'); 
Potgs=result.res;
clear result
result = load(sprintf(POT_FILE_FORMAT_SPHERE,'Potdm', MPc),'res');
Potdm = result.res;
clear result
result = load(sprintf(POT_FILE_FORMAT_SPHERE,'Potdot', MPc),'res');
Potdot = result.res;

