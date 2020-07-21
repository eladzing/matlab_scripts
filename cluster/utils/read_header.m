function [aexp h om0 ol0 ob0 ng dx] = read_header(MPSec)
% reads the file header for the box parameters. 
% useful for finding the cosmological parameters

% @param MPSec The required Mpc resolution (1, 2, 4, or 8)
% @returns  A 256^3 matrix of singles.
%

narginchk(1,1);

global FILE_FORMAT;
result = read_cube_header(sprintf(FILE_FORMAT, 'rhodm', MPSec));
aexp = result.aexpn;
h = result.hubble;
om0 = result.om0;
ol0 = result.ol0;
ob0 = result.ob0;
ng = result.nx;
dx = result.dx;