function res= redshift2snap(zred )
%REDSHIFT2SNAP get the snapshot closest to the given redshift
%   Loads pre-made table of snapshot to redshift and returns the right
%   redshift


global snapRedshifts

res=zeros(size(zred));

for i=1:length(zred)
    
    [~,ind]=min(abs(snapRedshifts-zred(i)));
    
    res(i)=ind-1;
end

end

