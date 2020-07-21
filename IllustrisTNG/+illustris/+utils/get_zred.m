function zred = get_zred(snap,simPath)
%GET_ZRED - get the redshift for a given snapshot
%   load the header fo the first chunk and extract the redshift, bp is the
%   optional basepath of the simulation

global BASEPATH

pat=BASEPATH;

if exist('simPath','var')
    pat=simPath;
end

header=illustris.groupcat.loadHeader(pat,snap,0);

zred=header.Redshift;
if zred<1e-10
    zred=0;
end

end

