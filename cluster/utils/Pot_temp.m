function result = Pot_temp(name)
% Utility function to load the potential cube
%
% @param MPSec The required Mpc resolution (1, 2, 4, or 8)
% @returns  A 256^3 matrix of singles.
%

error(nargchk(1,1,nargin,'struct'));

%global FILE_FORMAT;
ff='/home/eladzing/work/sshfs/sungate1/poteng/output/%s.dat';
%ff='/home/eladzing/work/clusters/datacubes/energy/potential/output/%s.dat';
result = read_cube(sprintf(ff, name));
result = result.data;
