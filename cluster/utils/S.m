function result = S(MPSec)
% Utility function to calculate the entropy cube according to:
% S = T/(\rho_{gas}^{2/3})
%
% @param MPSec The required Mpc resolution (1, 2, 4, or 8)
% @returns  A 256^3 matrix of singles.
%

narginchk(1,1);

result = T(MPSec)./((RHOG(MPSec)).^(2/3)); % in units K/(Msun^(2/3)/Mpc^2 

% use f_ent form units.m to change to keV cm^2 
