function [pgs,pdm,pdot] = Pot_temp_all(box)
% Utility function to load the potential cube
%
% @param MPSec The required Mpc resolution (1, 2, 4, or 8)
% @returns  A 256^3 matrix of singles.
%

error(nargchk(1,1,nargin,'struct'));

global POT_FILE_FORMAT;
%ff='/home/eladzing/work/sshfs/sungate1/poteng/output/%s.dat';
%ff='/home/eladzing/work/clusters/datacubes/energy/potential/data/%s/.dat';
result = read_cube(sprintf(POT_FILE_FORMAT,'potgas',box));
pgs = result.data;

fir j=result = read_cube(sprintf(POT_FILE_FORMAT,'potdm',box));
pdm = result.data;

result = read_cube(sprintf(POT_FILE_FORMAT,'potdot',box));
pdot = result.data;

