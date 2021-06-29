function result = MAKE_MASS_WEIGHTED_VAR_PROFILE(data, rho_sp)
% Utility function to generate a mass-weighted variance profile from 
% a given spherical cube (See: cart2sphere)
% Same as: MAKE_MASS_WEIGHTED_PROFILE(data.^2) - MAKE_MASS_WEIGHTED_PROFILE(data).^2
%
% @param data    The spherical data cube.
% @param rho_sp  The density profile in a spherical data cube.
%
% @returns     The mass-weighted variance profile.
%
error(nargchk(2,2,nargin,'struct'));

result = MAKE_MASS_WEIGHTED_PROFILE(data.^2, rho_sp) - MAKE_MASS_WEIGHTED_PROFILE(data, rho_sp).^2;