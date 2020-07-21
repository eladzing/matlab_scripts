function result = MAKE_PROFILE(data)
% Utility function to generate a volume-weighted profile from 
% a given spherical cube (See: cart2sphere)
%
% @param data  The spherical data cube.
%
% @returns     The volume-weighted profile.
%
error(nargchk(1,1,nargin,'struct'));

result = squeeze(sum(sum(data,3),2))'/(size(data,2) * size(data,3));