function result = RHODM(MPSec)
% Utility function to load the dark matter density cube
%
% @param MPSec The required Mpc resolution (1, 2, 4, or 8)
% @returns  A 256^3 matrix of singles.
%

error(nargchk(1,1,nargin,'struct'));

global FILE_FORMAT;
result = read_cube(sprintf(FILE_FORMAT, 'rhodm', MPSec));
result = result.data;