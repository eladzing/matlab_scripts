function [pdm pgs pdot ro]=poteng_temp(name,box)

%error(nargchk(1,1,nargin,'struct'));

%global FILE_FORMAT;
%ff='/home/eladzing/work/sshfs/sungate1/poteng/output/%s.dat';
%ff='/home/eladzing/work/clusters/datacubes/energy/potential/output/%s.dat';
global FILE_FORMAT


fname=sprintf('%sdm',name)

result = read_cube(sprintf(FILE_FORMAT,fname,box));
pdm = result.data;

fname=sprintf('%sgas',name)

result = read_cube(sprintf(FILE_FORMAT,fname,box));
pgs = result.data;

fname=sprintf('%sdot',name)

result = read_cube(sprintf(FILE_FORMAT,fname,box));
pdot = result.data;

fname='rhog';
result = read_cube(sprintf(FILE_FORMAT,fname,box));
ro = result.data;

