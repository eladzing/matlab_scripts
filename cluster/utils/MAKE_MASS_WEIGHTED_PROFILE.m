function result = MAKE_MASS_WEIGHTED_PROFILE(data, rho_sp)
% Utility function to generate a mass-weighted profile from 
% a given spherical cube (See: cart2sphere)
%
% @param data    The spherical data cube.
% @param rho_sp  The density profile in a spherical data cube.
%
% @returns     The mass-weighted profile.
%
error(nargchk(2,2,nargin,'struct'));

result = squeeze(sum(sum(data.*rho_sp,3),2))'./squeeze(sum(sum(rho_sp,3),2))';