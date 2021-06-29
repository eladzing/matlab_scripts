function result = read_galCube(filename)
% Utility function to load a uniform cube of galaxy data mde by Nir  
%
% @param filename  The file to load.
% @returns  A 300^3 matrix of singles.


  narginchk(1,1);

fid = fopen(filename, 'rb','ieee-be');
ntemp = fread(fid,1, 'int');
Ngrid  = fread(fid, 1, 'int');
Lgrid  = fread(fid, 1, 'float32');
Grid   = fread(fid, Ngrid*Ngrid*Ngrid, 'float32');

fclose(fid);

global NCELL
global boxSize
global RVIR

NCELL=Ngrid;
boxSize = 2*Lgrid;
RVIR = Lgrid/2;

result.Ngrid=Ngrid;
result.Lgrid=Lgrid;
result.Grid=reshape(Grid,[Ngrid Ngrid Ngrid]);
end
