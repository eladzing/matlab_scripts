function gasStruct =  addTemperature( gasStruct )
%ADDTEMP calculate temperature and add as array to gas structure
%   Makes use of the calcTemperature function 

    if isfield(gasStruct,'ElectronAbundance')
        gasStruct.Temperature = illustris.utils.calcTemperature(double(gasStruct.InternalEnergy),double(gasStruct.ElectronAbundance));
    else
        
        gasStruct.Temperature = illustris.utils.calcTemperature(double(gasStruct.InternalEnergy));
    end


end

