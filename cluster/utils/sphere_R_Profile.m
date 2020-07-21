function result = sphere_R_Profile(MPSec)
% Generates the 'canonical' radius profile for the spherical coordinates
% system.
%
% @param MPSec  The required Mpc resolution (1, 2, 4, or 8)
%
% @returns  The radius profile
%
error(nargchk(1,1,nargin,'struct'));

h = 0.7;
result = (0.5:0.5:128)*(MPSec/h)/256;