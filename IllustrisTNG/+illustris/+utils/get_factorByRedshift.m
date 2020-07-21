function res = get_factorByRedshift(name,zred)
%GET_REDSHIFTFACTOR get the redshift dependant unit factor of the Illustris simmulation units.
%   Detailed explanation goes here

global illUnits

snap0=illUnits.snap;


snaps=illustris.utils.redshift2snap(zred);
res=zeros(size(zred));

for i=1:length(snaps)

    illustris.utils.set_illUnits(snaps(i));
    
    res(i)=illUnits.(name);
end

illustris.utils.set_illUnits(snap0);

end

