function res= snap2aexpn(snap )
%SNAP2AEXPN get the expansion parameter of a given snapshot number
%   Loads pre-made table of snapshot to redshift and returns the right
%   redshift

global snapRedshifts

zred=snapRedshifts(snap+1);

res=1./(zred+1);

end

