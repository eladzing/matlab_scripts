function [Ek Ekr] = Ekin(MPc)
% Utility function to load the kinetic energy cube cube i
%
% @param MPSec The required Mpc resolution (1, 2, 4, or 8)
% @returns  A 256^3 matrix of singles.
%

error(nargchk(1,1,nargin,'struct'));

global FILE_FORMAT_MAT
result = load(sprintf(FILE_FORMAT_MAT,'Ekin',MPc),'res'); 
result = single(result.res);
Ek=result(:,:,:,1); 
Ekr=result(:,:,:,2);
