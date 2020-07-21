function result = MAKE_VAR_PROFILE(data)
% Utility function to generate a volume-weighted variance profile from 
% a given spherical cube (See: cart2sphere)
% Same as: MAKE_PROFILE(data.^2) - MAKE_PROFILE(data).^2
%
% @param data  The spherical data cube.
%
% @returns     The volume-weighted variance profile.
%
error(nargchk(1,1,nargin,'struct'));

result = MAKE_PROFILE(data.^2) - MAKE_PROFILE(data).^2;
