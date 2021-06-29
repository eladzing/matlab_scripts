function result = Pot(MPSec)
% Utility function to load the potential cube
%
% @param MPSec The required Mpc resolution (1, 2, 4, or 8)
% @returns  A 256^3 matrix of singles.
%

narginchk(1,1);

%global FILE_FORMAT;
ff='/home/eladzing/work/sshfs/sungate1/poteng/data/CL6/CSF/%s_a1.001L%dMpc.dat';
result = read_cube(sprintf(ff, 'pot', MPSec));
result = result.data;
