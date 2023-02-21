function gasStruct =  addCellRadius( gasStruct,unitType )
%ADDCEELLRADIUS calculate the 'spherical' cell radius returns in adjusted
%units
%   assume cell is a sphere and find the radius based on mass and density.
%    mass and density are in simulation units (comoving /h) by default but 
%    can be set to physical units (no /h) .

if ~exist('unitType','var')
    unitType='simulation'; % comoving /h - same as in the structure
end
global illUnits

switch lower(unitType)
    case 'simulation'
        gasStruct.CellRadius= (3/(4*pi).* ( gasStruct.Masses)./ (gasStruct.Density)).^(1/3);
    case 'physical'
        gasStruct.CellRadius= (3/(4*pi).* ( gasStruct.Masses.*illUnits.massUnit)./ (gasStruct.Density.*illUnits.densityUnit)).^(1/3);
    otherwise
        error('%s - unknown unit type: %s',current_function().upper,unitType);
end


end

