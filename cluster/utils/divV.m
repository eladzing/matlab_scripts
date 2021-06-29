function div = divV(MPc)
% Utility function to load the velocity divergence field
% The units are a bit funny: [u]=km/sec [div]=1/Mpc proper 

% @param MPSec The required Mpc resolution (1, 2, 4, or 8)
% @returns  A 256^3 matrix of singles.
%

error(nargchk(1,1,nargin,'struct'));

global FILE_FORMAT_MAT
result = load(sprintf(FILE_FORMAT_MAT,'divV',MPc),'res'); 
result = single(result.res);
div=result;

