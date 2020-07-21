function result = S_sphere(MPSec)
% Utility function to load the entropy cube in spherical coordinates
% (See: cart2sphere)
%
% @param MPSec The required Mpc resolution (1, 2, 4, or 8)
% @returns  A 256^3 matrix of singles.
%

narginchk(1,1);

result = T_sphere(MPSec)./((RHOG_sphere(MPSec)).^(2/3));

%global FILE_FORMAT_SPHERE;
%result = load(sprintf(FILE_FORMAT_SPHERE,'S', MPc),'res'); 
%result = single(result.res);
