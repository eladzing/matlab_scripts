function res = get_stellar_mass(subs,mask)
%GET_STELLAR_MASS - get the galaxy stellar mass for TNG gal's, i.e. 
%    stellar mass within 2*rhalf. 

if ~exist('mask','var')
    mask=true(1,subs.count);
end

global illUnits
res = double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,mask).*illUnits.massUnit); 

end

