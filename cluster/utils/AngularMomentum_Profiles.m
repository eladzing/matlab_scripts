function [Jx_Profile Jy_Profile Jz_Profile J_Profile Jtheta_Profile Jphi_Profile] = AngularMomentum_Profiles(Jx,Jy,Jz)
% Calculates the angular momentum profiles from the angular momentum cubes.
%
% @param Jx,Jy,Jz  The angular momentum cubes in the X, Y, and Z directions
%                  respectively.
% @returns  The angular momentum profiles in cartesian and spherical
%           coordinates.
%
error(nargchk(3,3,nargin,'struct'));

% Cartesian coordinates
Jx_Profile = MAKE_PROFILE_FROM_CUBE(Jx);
Jy_Profile = MAKE_PROFILE_FROM_CUBE(Jy);
Jz_Profile = MAKE_PROFILE_FROM_CUBE(Jz);

% Spherical coordinates
J_Profile      = sqrt(Jx_Profile.^2 + Jy_Profile.^2 + Jz_Profile.^2);
Jtheta_Profile = acos(Jz_Profile ./ J_Profile);
Jphi_Profile   = atan2(Jy_Profile,Jx_Profile);