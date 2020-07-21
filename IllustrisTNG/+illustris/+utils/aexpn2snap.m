function res= aexpn2snap(aexpn )
%REDSHIFT2SNAP get the snapshot closest to the given redshift
%   Loads pre-made table of snapshot to redshift and returns the right
%   redshift

zred=1./aexpn-1;

global DEFAULT_MATFILE_DIR
% load mat file
load([DEFAULT_MATFILE_DIR '/snaps2redshift.mat']);

for i=1:length(zred)
    
    [~,ind]=min(abs(snaps2redshift-zred(i)));
    
    res(i)=ind-1;
end

end

