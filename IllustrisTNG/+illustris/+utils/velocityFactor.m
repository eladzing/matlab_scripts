function unitFac = velocityFactor(snap,velType)
%VELOCITYFACTOR get the right unit factor for the different velocities in
%Illustrsis
%   different velocities have different units in the Illustris catalogs -
%   FoF, subfind, particles etc. have different units in velocity. This
%   function supplies the right conversion factor

if snap<=0
    error('%s - snap must be positive: %s',current_function().upper,snap);
end

if floor(snap)~=snap
    error('%s - snap should be an integer: %s',current_function().upper,snap);
end


switch(lower(velType))
    case{'gas','dm','darkmatter','particle','part','particles','bh','stars'}
        unitFac=sqrt(illustris.utils.snap2aexpn(snap));
    case{'host','halo','halos','fof'}
            unitFac=1./(illustris.utils.snap2aexpn(snap));
    case{'subfind','sub','subs','subhalo','subhalos'}
    unitFac=1;
    otherwise
        error('velocityFactor -Illegal argument for velType: %s',velType);
end



end

