function result = RHOTOT(MPSec)
% Utility function to load the total density cube
%
% @param MPSec The required Mpc resolution (1, 2, 4, or 8)
% @returns  A 256^3 matrix of singles.
%

narginchk(1,1);

global FILE_FORMAT;
result = read_cube(sprintf(FILE_FORMAT, 'rhotot', MPSec));
result = result.data;