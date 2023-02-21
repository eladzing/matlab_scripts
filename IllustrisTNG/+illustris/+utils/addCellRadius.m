function gasStruct =  addCellRadius( gasStruct )
%ADDCEELLRADIUS calculate the 'spherical' cell radius returns in adjusted
%units
%   assume cell is a sphere and find the radius based on mass and density. 
%    mass and density are in adjusted (no /h)  physical units. 

global illUnits

gasStruct.CellRadius= (3/(4*pi).* ( gasStruct.Masses.*illUnits.massUnit)./ (gasStruct.Density.*illUnits.densityUnit)).^(1/3);

end

