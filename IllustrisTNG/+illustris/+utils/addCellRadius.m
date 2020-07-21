function gasStruct =  addCellRadius( gasStruct )
%ADDTEMP calculate temperature and add as array to gas structure
%   Makes use of the calcTemperature function

global illUnits

gasStruct.CellRadius= (3/(4*pi).* ( gasStruct.Masses.*illUnits.massUnit)./ (gasStruct.Density.*illUnits.densityUnit)).^(1/3);

end

