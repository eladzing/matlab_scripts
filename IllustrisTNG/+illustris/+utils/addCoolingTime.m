function gasStruct =  addCoolingTime( gasStruct )
%ADDTEMP calculate temperature and add as array to gas structure
%   Makes use of the calcTemperature function 

 gasStruct.CoolingTime = illustris.utils.calcCoolingTime(double(gasStruct.Density),...
     double(gasStruct.InternalEnergy),double(gasStruct.GFM_CoolingRate));


end

