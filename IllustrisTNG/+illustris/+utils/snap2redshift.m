function res= snap2redshift(snap )
%SNAP2REDSHIFT get the redshift of a given snapshot number
%   Loads pre-made table of snapshot to redshift and returns the right
%   redshift

global snapRedshifts

res=snapRedshifts(snap+1);


end

