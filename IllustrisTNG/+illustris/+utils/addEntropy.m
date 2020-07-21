function gasStruct =  addEntropy( gasStruct )
%ADDEntropy calculate entropy and add as array to gas structure
%   Makes use of the calcEntropy function. if Temperature ahas not yet been
%   defined then it will be....

if ~isfield(gasStruct,'Temperature')
    gasStruct=illustris.utils.addTemperature(gasStruct);
end

gasStruct.Entropy = illustris.utils.calcEntropy(gasStruct.Temperature,double(gasStruct.Density));
end

