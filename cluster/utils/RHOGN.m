function result = RHOGN(MPSec)
% Utility function to load the gas density cube
%
% @param MPSec The required Mpc resolution (1, 2, 4, or 8)
% @returns  A 256^3 matrix of singles.
%
units;
fac=Ms/mm/Mpc^3;
narginchk(1,1);

global FILE_FORMAT;
result = read_cube(sprintf(FILE_FORMAT, 'rhog', MPSec));
result = result.data.*fac;