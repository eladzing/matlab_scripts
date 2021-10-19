function res = get_stellar_mass(subs,varargin) %mask,type)
%GET_STELLAR_MASS - get the galaxy stellar mass for TNG gal's, i.e.
%    stellar mass within 2*rhalf.

type='gal'; % stellar mass within galactic radius
mask=true(1,subs.count);

i=1;
while i<=length(varargin)
    switch lower(varargin{i})
        case 'mask'
            i=i+1;
            mask=varargin{i};
        case {'full','sub'}
            type='full';  % all stellar mass in the subhalo
        case {'gal','galaxy'}
            type='gal';
        case 'cgm'  % all stellar mass in the subhalo, beyond the galactic radius
            type='cgm';
    end
    i=i+1;
end


global illUnits
switch type
    case 'gal'
        res = double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,mask).*illUnits.massUnit);
    case 'full'
        res = double(subs.SubhaloMassType(illustris.partTypeNum('stars')+1,mask).*illUnits.massUnit);
    case 'cgm='
        res = (double(subs.SubhaloMassType(illustris.partTypeNum('stars')+1,mask))-...
            double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,mask))).*illUnits.massUnit;
        
end

