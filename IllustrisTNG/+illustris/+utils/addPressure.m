function gasStruct =  addPressure( gasStruct )
%ADDTEMP calculate temperature and add as array to gas structure
%   Makes use of the calcTemperature function 

global illUnits

if isfield(gasStruct,'Temperature')
    tmp=gasStruct.Temperature;
else
    
    tmp=illustris.utils.calcTemperature(double(gasStruct.InternalEnergy),double(gasStruct.ElectronAbundance));
    
end

% convert density to number density 
nDensity=double(gasStruct.Density.*illUnits.numberDensityFactor); 

gasStruct.Pressure = illustris.utils.calcPressure(nDensity,tmp);

end

